// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_FUNCTIONS_INTERPOLATE_HH
#define DUNE_DPG_FUNCTIONS_INTERPOLATE_HH

#include <algorithm>
#include <dune/functions/common/functionfromcallable.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/geometry/referenceelements.hh>

namespace Dune {

template<typename CoarseBasis, typename FineBasis, typename VectorType>
VectorType interpolateToUniformlyRefinedGrid(
    const CoarseBasis& coarseBasis,
    const FineBasis& fineBasis,
    const VectorType& v_coarse)
{
  auto coarseGridView = coarseBasis.gridView();

  auto coarseLocalView = coarseBasis.localView();
  auto fineLocalView = fineBasis.localView();

  auto coarseLocalIndexSet = coarseBasis.localIndexSet();
  auto fineLocalIndexSet = fineBasis.localIndexSet();

  VectorType v_fine(fineBasis.size()); v_fine = 0;

  for(const auto& e : elements(coarseGridView)) {
    coarseLocalView.bind(e);
    coarseLocalIndexSet.bind(coarseLocalView);

    VectorType local_v_coarse(coarseLocalView.size());

    for (size_t i=0, nCoarse=coarseLocalView.size(); i<nCoarse; i++)
    {
      auto row = coarseLocalIndexSet.index(i)[0];
      local_v_coarse[i] = v_coarse[row];
    }

    // TODO: Seems, I have to reed the Dune documentation some more
    auto& coarseLocalBasis = coarseLocalView.tree().finiteElement().localBasis();
    using CoarseFiniteElement = std::decay_t<decltype(coarseLocalView.tree().finiteElement())>;
    using CoarseLocalBasis = typename CoarseFiniteElement::Traits::LocalBasisType;
    using FiniteElementRange = typename CoarseLocalBasis::Traits::RangeType;
    using FunctionBaseClass = typename Dune::LocalFiniteElementFunctionBase<CoarseFiniteElement>::type;
    using LocalDomain = typename CoarseLocalBasis::Traits::DomainType;
    constexpr unsigned int dim = CoarseLocalBasis::Traits::dimDomain;

    auto globalGeometry = e.geometry();
    assert(globalGeometry.affine());

    for (const auto& subE : descendantElements(e, 1))
    {
      fineLocalView.bind(subE);
      fineLocalIndexSet.bind(fineLocalView);

      auto globalSubGeometry = subE.geometry();

      // TODO: replace double with correct typedef
      const auto& subGeometry = AffineGeometry<double, dim, dim>
          (subE.type(),
           globalGeometry.local(
             globalSubGeometry.global(ReferenceElements<double, dim>
               ::general(subE.type()).position(0,dim))),
           globalSubGeometry.jacobianTransposed({})
             .leftmultiply(globalGeometry.jacobianTransposed({})));
      auto localF = [subGeometry, &coarseLocalBasis, local_v_coarse]
                    (const LocalDomain& x) {
        auto xCoarse = subGeometry.global(x);

        std::vector<FiniteElementRange> shapeValues;
        coarseLocalBasis.evaluateFunction(xCoarse, shapeValues);
        FiniteElementRange y
          = std::inner_product(shapeValues.cbegin(), shapeValues.cend(),
                               local_v_coarse.begin(), 0.);
        return y;
      };

      using FunctionFromCallable = typename Dune::Functions::FunctionFromCallable<FiniteElementRange(LocalDomain), decltype(localF), FunctionBaseClass>;

      auto interpolationValues = std::vector<FiniteElementRange>();
      fineLocalView.tree().finiteElement().localInterpolation()
          .interpolate(FunctionFromCallable(localF), interpolationValues);

      for (size_t i=0, nFine=fineLocalView.size(); i<nFine; i++)
      {
        auto row = fineLocalIndexSet.index(i)[0];
        v_fine[row] = interpolationValues[i];
      }
    }
  }

  return v_fine;
}

} // end namespace Dune

#endif
