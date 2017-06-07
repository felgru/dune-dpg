// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_APPROXIMATE_SCATTERING_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_APPROXIMATE_SCATTERING_HH

#include <memory>
#include <tuple>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/localevaluation.hh>
#include <dune/dpg/quadrature.hh>

#include <dune/istl/bvector.hh>

#include "svdkernelapproximation.hh"
#include "waveletkernelapproximation.hh"

namespace Dune {

/**
 * \brief This constructs the right hand side vector of a DPG system.
 *
 * \tparam SolutionSpaces  tuple of solution spaces
 */
template<class SolutionSpace,
         class KernelApproximation>
class ApproximateScatteringAssembler
{
public:
  enum : unsigned int { dim = SolutionSpace::GridView::dimension };
  using Direction = FieldVector<double, dim>;

  ApproximateScatteringAssembler () = delete;
  /**
   * \brief constructor for ApproximateScatteringAssembler
   *
   * \param solutionSpace    a reference to a solution space (on the host grid)
   * \param kernel           an approximation of the scattering kernel
   *
   * \note For your convenience, use make_ApproximateScatteringAssembler()
   *       instead.
   */
  ApproximateScatteringAssembler (
      const SolutionSpace& solutionSpace,
      const KernelApproximation& kernel)
             : solutionSpace(solutionSpace),
               kernelApproximation(kernel)
  {}

  /**
   * \brief Assemble a vector representing the scattering integral
   *        in the solution space basis for a given set of
   *        discrete ordinate solutions.
   *
   * \todo For now, the scattering kernel is assumed to be constant
   *       in space and normed to integral 1.
   *
   * \param[out] scattering  the scattering vector
   * \param[in]  x           the vectors of the solutions of the
   *                           previous iteration
   * \param[in]  si          index of scattering direction (0 < si < numS)
   */
  void precomputeScattering
          (BlockVector<FieldVector<double,1> >& scattering,
           const std::vector<BlockVector<FieldVector<double,1> >>& x,
           size_t si);

private:
  const SolutionSpace& solutionSpace;
  const KernelApproximation& kernelApproximation;
};

/**
 * \brief Creates a ScatteringAssembler for a DPG discretization,
 *        deducing the target type from the types of arguments.
 *
 * \param  solutionSpace    a reference to a solution spaces (on host grid).
 *                          It is assumed that the solutions given later have
 *                          all been interpolated to this solutionSpace.
 * \param  kernel           an approximation of the scattering kernel
 */
template<class SolutionSpace,
         class KernelApproximation>
auto make_ApproximateScatteringAssembler(
      const SolutionSpace& solutionSpace,
      const KernelApproximation& kernel)
     -> ApproximateScatteringAssembler<SolutionSpace, KernelApproximation>
{
  return ApproximateScatteringAssembler<SolutionSpace, KernelApproximation>
                            (solutionSpace, kernel);
}

namespace detail {

template <class SolutionSpace>
struct GetLocalApproximateScattering
{
using SolutionLocalView = typename SolutionSpace::LocalView;
using SolutionLocalIndexSet = typename SolutionSpace::LocalIndexSet;

template <class VectorType,
          class Element,
          class KernelApproximation>
inline static void interiorImpl(
    const SolutionLocalView& solutionLocalView,
    VectorType& localScattering,
    VectorType& localNormSquared,
    const SolutionLocalIndexSet& solutionLocalIndexSet,
    // unsigned int quadratureOrder,
    const Element& element,
    const KernelApproximation& kernelApproximation,
    const std::vector<BlockVector<FieldVector<double,1> >>& x,
    size_t si)
{
  const int dim = Element::mydimension;
  auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& solutionLocalFiniteElement
      = solutionLocalView.tree().finiteElement();

  assert(localScattering.size() == solutionLocalFiniteElement.size());
  assert(localNormSquared.size() == solutionLocalFiniteElement.size());

  /* TODO:
   * - Adapt quadrature also to the kernel k
   */
  const unsigned int quadratureOrder
      = 2 * solutionLocalFiniteElement.localBasis().order();

  typename detail::ChooseQuadrature<SolutionSpace, SolutionSpace, Element>::type quad
    = detail::ChooseQuadrature<SolutionSpace, SolutionSpace, Element>
      ::Quadrature(element, quadratureOrder, nullptr);

  // Loop over all quadrature points
  for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPos);

    std::vector<FieldVector<double,1>> shapeFunctionValues;
    solutionLocalFiniteElement.localBasis().
        evaluateFunction(quadPos, shapeFunctionValues);

    const size_t numS = x.size();
    Eigen::VectorXd uValues(numS);

    for (size_t scatteringAngle=0;
         scatteringAngle<numS; ++scatteringAngle) {
      double uValue = 0; // in direction of scatteringAngle
      // Evaluate all shape function values at this point
      std::vector<FieldVector<double,1>> shapeFunctionValues;
      solutionLocalFiniteElement.localBasis().
          evaluateFunction(quadPos, shapeFunctionValues);
      for (size_t j=0, jMax=shapeFunctionValues.size(); j<jMax; j++)
      {
        /* This assumes that solutionLocalIndexSet and
         * doesn't change for different * scattering angles.
         */
        auto row =
            solutionLocalIndexSet.index(j)[0];
        uValue += x[scatteringAngle][row] * shapeFunctionValues[j];
      }
      uValues(scatteringAngle) = uValue;
    }

    kernelApproximation.applyToVector(uValues);
    /* TODO: shouldn't we integrate over the angles somewhere? */

    const double integrationWeight = quad[pt].weight() * integrationElement;
    const double factor = uValues(si) * integrationWeight;
    for (size_t i=0, i_max=localScattering.size(); i<i_max; i++) {
      localScattering[i] += factor * shapeFunctionValues[i];
      localNormSquared[i] += shapeFunctionValues[i] * shapeFunctionValues[i]
                           * integrationWeight;
    }
  }
}

};

} // end namespace detail


template<class SolutionSpace,
         class KernelApproximation>
void ApproximateScatteringAssembler<SolutionSpace, KernelApproximation>::
precomputeScattering(BlockVector<FieldVector<double,1> >& scattering,
                     const std::vector<BlockVector<FieldVector<double,1> >>& x,
                     size_t si)
{
  using namespace Dune::detail;

  typedef typename SolutionSpace::GridView GridView;
  GridView gridView = solutionSpace.gridView();

  scattering.resize(solutionSpace.size());
  scattering = 0;
  BlockVector<FieldVector<double,1>> normSquared(scattering.size());
  normSquared = 0;

  // Views on the FE bases on a single element
  auto solutionLocalView = solutionSpace.localView();
  auto solutionLocalIndexSet = solutionSpace.localIndexSet();

  for(const auto& e : elements(gridView)) {

    // Bind the local FE basis view to the current element
    solutionLocalView.bind(e);
    solutionLocalIndexSet.bind(solutionLocalView);

    // Now get the local contribution to the scattering functional

    BlockVector<FieldVector<double,1>>
        localScattering(solutionLocalView.size());
    localScattering = 0;
    BlockVector<FieldVector<double,1>> localNormSquared(localScattering.size());
    localNormSquared = 0;

    detail::GetLocalApproximateScattering<SolutionSpace>::interiorImpl(
            solutionLocalView,
            localScattering,
            localNormSquared,
            solutionLocalIndexSet,
            // unsigned int quadratureOrder,
            e,
            kernelApproximation,
            x,
            si);

    using MultiIndex
        = typename std::decay_t<decltype(solutionLocalIndexSet)>::MultiIndex;
    iterateOverLocalIndexSet(
        solutionLocalIndexSet,
        [&](size_t i, MultiIndex gi)
        {
          scattering[gi[0]]
              += localScattering[i];
        },
        [](size_t i){},
        [&](size_t i, MultiIndex gi, double wi)
        {
          scattering[gi[0]]
              += wi * localScattering[i];
        }
    );

    iterateOverLocalIndexSet(
        solutionLocalIndexSet,
        [&](size_t i, MultiIndex gi)
        {
          normSquared[gi[0]]
              += localNormSquared[i];
        },
        [](size_t i){},
        [&](size_t i, MultiIndex gi, double wi)
        {
          normSquared[gi[0]]
              += wi * localNormSquared[i];
        }
    );
  }

  for(size_t i=0, iMax=scattering.size(); i<iMax; ++i)
    scattering[i] /= normSquared[i];
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_APPROXIMATE_SCATTERING_HH
