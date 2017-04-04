// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_APPROXIMATE_SCATTERING_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_APPROXIMATE_SCATTERING_HH

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
template<class SolutionSpaces,
         class KernelApproximation>
class ApproximateScatteringAssembler
{
public:
  enum : unsigned int {
    dim = std::tuple_element_t<0, SolutionSpaces>::GridView::dimension };
  using Direction = FieldVector<double, dim>;

  ApproximateScatteringAssembler () = delete;
  /**
   * \brief constructor for ApproximateScatteringAssembler
   *
   * \param solutionSpaces   a tuple of solution spaces
   * \param kernel           an approximation of the scattering kernel
   * \param si  index of the scattering direction
   *
   * \note For your convenience, use make_ApproximateScatteringAssembler()
   *       instead.
   */
  ApproximateScatteringAssembler (const SolutionSpaces& solutionSpaces,
                                  const KernelApproximation& kernel,
                                  size_t si)
             : solutionSpaces(solutionSpaces),
               si(si),
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
   */
  template<size_t solutionSpaceIndex>
  void precomputeScattering
          (BlockVector<FieldVector<double,1> >& scattering,
           const std::vector<BlockVector<FieldVector<double,1> >>& x);

private:
  SolutionSpaces solutionSpaces;
  const size_t si;
  const KernelApproximation& kernelApproximation;
};

/**
 * \brief Creates a ScatteringAssembler for a DPG discretization,
 *        deducing the target type from the types of arguments.
 *
 * \param  solutionSpaces   a tuple of solution spaces
 * \param  kernel           an approximation of the scattering kernel
 * \param  si               index of scattering direction (0 < si < numS)
 */
template<class SolutionSpaces,
         class KernelApproximation>
auto make_ApproximateScatteringAssembler(
      const SolutionSpaces& solutionSpaces,
      const KernelApproximation& kernel,
      size_t si)
     -> ApproximateScatteringAssembler<SolutionSpaces, KernelApproximation>
{
  return ApproximateScatteringAssembler<SolutionSpaces, KernelApproximation>
                            (solutionSpaces, kernel, si);
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
    size_t globalSolutionSpaceOffset,
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
        /* This assumes that solutionLocalIndexSets and
         * globalSolutionSpaceOffset don't change for different
         * scattering angles.
         */
        auto row =
            solutionLocalIndexSet.index(j)[0]
          + globalSolutionSpaceOffset;
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


template<class SolutionSpaces,
         class KernelApproximation>
template<size_t solutionSpaceIndex>
void ApproximateScatteringAssembler<SolutionSpaces, KernelApproximation>::
precomputeScattering(BlockVector<FieldVector<double,1> >& scattering,
                     const std::vector<BlockVector<FieldVector<double,1> >>& x)
{
  using namespace Dune::detail;

  typedef typename std::tuple_element<0,SolutionSpaces>::type::GridView
          GridView;
  GridView gridView = std::get<0>(solutionSpaces).gridView();

  /* set up global offsets */
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

  computeOffsets(globalSolutionSpaceOffsets, solutionSpaces);

  const size_t globalSolutionSpaceOffset =
      globalSolutionSpaceOffsets[solutionSpaceIndex];

  scattering.resize(std::get<solutionSpaceIndex>(solutionSpaces).size());
  scattering = 0;
  BlockVector<FieldVector<double,1>> normSquared(scattering.size());
  normSquared = 0;

  // Views on the FE bases on a single element
  auto solutionLocalViews = getLocalViews(solutionSpaces);
  auto solutionLocalIndexSets = getLocalIndexSets(solutionSpaces);

  for(const auto& e : elements(gridView)) {

    // Bind the local FE basis view to the current element
    /* TODO: only bind the space we use later */
    bindLocalViews(solutionLocalViews, e);
    bindLocalIndexSets(solutionLocalIndexSets, solutionLocalViews);

    // Now get the local contribution to the scattering functional

    BlockVector<FieldVector<double,1>>
        localScattering(std::get<solutionSpaceIndex>(solutionLocalViews)
                                                                  .size());
    localScattering = 0;
    BlockVector<FieldVector<double,1>> localNormSquared(localScattering.size());
    localNormSquared = 0;

    detail::GetLocalApproximateScattering
      < std::tuple_element_t<solutionSpaceIndex, SolutionSpaces>
      >::interiorImpl(
            std::get<solutionSpaceIndex>(solutionLocalViews),
            localScattering,
            localNormSquared,
            std::get<solutionSpaceIndex>(solutionLocalIndexSets),
            globalSolutionSpaceOffset,
            // unsigned int quadratureOrder,
            e,
            kernelApproximation,
            x,
            si);

    typename std::remove_reference<decltype(globalSolutionSpaceOffsets)>
                         ::type zeroes;
    for(size_t i = 0; i < std::tuple_size<SolutionSpaces>::value; i++)
      zeroes[i] = 0;
    auto scatteringCopier
        = localToGlobalRHSCopier
            <typename std::remove_reference<decltype(localScattering)>::type,
             typename std::remove_reference<decltype(scattering)>::type,
             typename std::remove_reference<decltype(solutionLocalIndexSets)>
                         ::type,
             typename std::remove_reference<decltype(zeroes)>::type>
            (localScattering, scattering,
             solutionLocalIndexSets,
             zeroes,
             zeroes);
    scatteringCopier(std::integral_constant<size_t, solutionSpaceIndex>{});

    auto normCopier
        = localToGlobalRHSCopier
            <typename std::remove_reference<decltype(localNormSquared)>::type,
             typename std::remove_reference<decltype(normSquared)>::type,
             typename std::remove_reference<decltype(solutionLocalIndexSets)>
                         ::type,
             typename std::remove_reference<decltype(zeroes)>::type>
            (localNormSquared, normSquared,
             solutionLocalIndexSets,
             zeroes,
             zeroes);
    normCopier(std::integral_constant<size_t, solutionSpaceIndex>{});
  }

  for(size_t i=0, iMax=scattering.size(); i<iMax; ++i)
    scattering[i] /= normSquared[i];
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_APPROXIMATE_SCATTERING_HH
