// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_FUNCTIONS_INTERPOLATE_HH
#define DUNE_DPG_FUNCTIONS_INTERPOLATE_HH

#include <algorithm>
#include <dune/functions/common/functionfromcallable.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>

namespace Dune {

/**
 * This function interpolates to a uniform refinement of a grid
 *
 * \pre The GridView of \p fineBasis has to differ from \p coarseBasis
 *      by exactly one level of uniform refinement.
 *
 * \param coarseBasis the GlobalBasis of the FE function on the coarser level
 * \param fineBasis   the GlobalBasis of the FE function on the finer level
 * \param v_coarse    the coefficient vector of the FE function on the
 *                    coarser level.
 * \return  the coefficient vector resulting from interpolating v_coarse
 *          to fineBasis.
 */
template<typename CoarseBasis, typename FineBasis, typename VectorType>
VectorType interpolateToUniformlyRefinedGrid(
    const CoarseBasis& coarseBasis,
    const FineBasis& fineBasis,
    const VectorType& v_coarse)
{
  auto coarseGridView = coarseBasis.gridView();

  auto coarseLocalView = coarseBasis.localView();
  auto fineLocalView = fineBasis.localView();

  VectorType v_fine(fineBasis.size()); v_fine = 0;

  for(const auto& e : elements(coarseGridView)) {
    coarseLocalView.bind(e);

    VectorType local_v_coarse(coarseLocalView.size());

    for (size_t i=0, nCoarse=coarseLocalView.size(); i<nCoarse; i++)
    {
      auto row = coarseLocalView.index(i)[0];
      local_v_coarse[i] = v_coarse[row];
    }

    auto& coarseLocalBasis = coarseLocalView.tree().finiteElement().localBasis();
    using CoarseFiniteElement = std::decay_t<decltype(coarseLocalView.tree().finiteElement())>;
    using CoarseLocalBasis = typename CoarseFiniteElement::Traits::LocalBasisType;
    using FiniteElementRange = typename CoarseLocalBasis::Traits::RangeType;
    using FunctionBaseClass = typename Dune::LocalFiniteElementFunctionBase<CoarseFiniteElement>::type;
    using LocalDomain = typename CoarseLocalBasis::Traits::DomainType;
    constexpr unsigned int dim = CoarseLocalBasis::Traits::dimDomain;

    auto globalGeometry = e.geometry();
    assert(globalGeometry.affine());

    for (const auto& subE : descendantElements(e, e.level()+1))
    {
      fineLocalView.bind(subE);

      auto globalSubGeometry = subE.geometry();

      // TODO: replace double with correct typedef
      const AffineGeometry<double, dim, dim> subGeometry
          (subE.type(),
           globalGeometry.local(
             globalSubGeometry.global(referenceElement<double, dim>
               (subE.type()).position(0,dim))),
           globalSubGeometry.jacobianTransposed({})
             .leftmultiply(globalGeometry.jacobianTransposed({})));
      auto localF = [&subGeometry, &coarseLocalBasis, &local_v_coarse]
                    (const LocalDomain& x) {
        auto xCoarse = subGeometry.global(x);

        std::vector<FiniteElementRange> shapeValues;
        coarseLocalBasis.evaluateFunction(xCoarse, shapeValues);
        FiniteElementRange y
          = std::inner_product(shapeValues.cbegin(), shapeValues.cend(),
                               local_v_coarse.begin(), FiniteElementRange{0.});
        return y;
      };

      using FunctionFromCallable = typename Dune::Functions::FunctionFromCallable<FiniteElementRange(LocalDomain), decltype(localF), FunctionBaseClass>;

      auto interpolationValues = std::vector<FiniteElementRange>();
      fineLocalView.tree().finiteElement().localInterpolation()
          .interpolate(FunctionFromCallable(localF), interpolationValues);

      for (size_t i=0, nFine=fineLocalView.size(); i<nFine; i++)
      {
        const auto row = fineLocalView.index(i)[0];
        v_fine[row] = interpolationValues[i];
      }
    }
  }

  return v_fine;
}

} // end namespace Dune

#endif
