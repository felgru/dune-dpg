// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_SCATTERING_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_SCATTERING_HH

#include <tuple>
#include <functional>
#include <memory>
#include <type_traits>

#include <dune/common/fvector.hh>

#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/quadrature.hh>

#include <dune/istl/bvector.hh>

namespace Dune {

/**
 * \brief This constructs the right hand side vector of a DPG system.
 *
 * \tparam TestSpaces      tuple of test spaces
 * \tparam SolutionSpaces  tuple of solution spaces
 */
template<class TestSpaces,
         class SolutionSpaces>
class ScatteringAssembler
{
public:
  ScatteringAssembler () = delete;
  /**
   * \brief constructor for ScatteringAssembler
   *
   * \note For your convenience, use make_ScatteringAssembler() instead.
   */
  constexpr ScatteringAssembler (const TestSpaces& testSpaces,
                                 const SolutionSpaces& solutionSpaces)
             : testSpaces(testSpaces),
               solutionSpaces(solutionSpaces)
  { };

  /**
   * \brief Assemble the vector corresponding to the scattering integral
   * for a given set of discrete ordinate solutions.
   *
   * \todo For now, the scattering kernel is assumed to be constant
   *       and normed to integral 1.
   *
   * \param[out] scattering  the scattering vector
   * \param[in]  x           the vectors of the solutions of the
   *                           previous iteration
   * \param[in] sVector
   * \param[in] kernelS
   */
  template<size_t solutionSpaceIndex, class Direction, class Function>
  void assembleScattering
         (BlockVector<FieldVector<double,1> >& scattering,
          const std::vector<BlockVector<FieldVector<double,1> >>& x,
          const std::vector<Direction>& sVector,
          Function& kernelS);

private:
  template<class Function, class Geometry, class QuadratureRule, class Direction>
  void evaluateKernelElement
    (std::vector<std::vector<double>>&,
    Function&,
    const Geometry&,
    const QuadratureRule&,
    const std::vector<Direction>&);

  TestSpaces     testSpaces;
  SolutionSpaces solutionSpaces;
};

/**
 * \brief Creates a ScatteringAssembler for a DPG discretization,
 *        deducing the target type from the types of arguments.
 *
 * \param  testSpaces       a tuple of test spaces
 * \param  solutionSpaces   a tuple of solution spaces
 */
template<class TestSpaces,
         class SolutionSpaces>
auto make_DPG_ScatteringAssembler(const TestSpaces& testSpaces,
                                  const SolutionSpaces& solutionSpaces)
     -> ScatteringAssembler<TestSpaces, SolutionSpaces>
{
  return ScatteringAssembler<TestSpaces, SolutionSpaces>
                            (testSpaces, solutionSpaces);
}

/**
 * \brief Creates a ScatteringAssembler for a saddlepoint discretization,
 *        deducing the target type from the types of arguments.
 *
 * \param  testSpaces       a tuple of test spaces
 * \param  solutionSpaces   a tuple of solution spaces
 */
template<class TestSpaces,
         class SolutionSpaces>
auto make_Saddlepoint_ScatteringAssembler(const TestSpaces& testSpaces,
                                          const SolutionSpaces& solutionSpaces)
     -> ScatteringAssembler<TestSpaces, SolutionSpaces>
{
  return ScatteringAssembler<TestSpaces, SolutionSpaces>
                            (testSpaces, solutionSpaces);
}

template<class TestSpaces,
         class SolutionSpaces>
template<class Function, class Geometry, class QuadratureRule, class Direction>
void ScatteringAssembler<TestSpaces, SolutionSpaces>::
     evaluateKernelElement(
                          std::vector<std::vector<double>>& k,
                          Function& kernelS,
                          const Geometry& localGeometry,
                          const QuadratureRule& quad,
                          const std::vector<Direction>& sVector)
{
  k.resize(quad.size());

  for ( size_t iQuad=0; iQuad < quad.size(); iQuad++ ) {

      k[iQuad].resize(sVector.size());
      // Position of the current quadrature point in the reference element
      auto quadPos = quad[iQuad].position();
      // Position of the current quadrature point in the current element
      auto mapQuadPos = localGeometry.global(quadPos);  // we get the global coordinate of quadPos
      // loop over the directions
      for(size_t iS=0; iS < sVector.size(); iS++)
      {
        k[iQuad][iS] = std::get<0>(kernelS)(mapQuadPos,sVector[iS]);
      }
  }
}


template<class TestSpaces,
         class SolutionSpaces>
template<size_t solutionSpaceIndex, class Direction, class Function>
void ScatteringAssembler<TestSpaces, SolutionSpaces>::
assembleScattering(BlockVector<FieldVector<double,1> >& scattering,
                   const std::vector<BlockVector<FieldVector<double,1> >>& x,
                   const std::vector<Direction>& sVector,
                   Function& kernelS)
{
  using namespace Dune::detail;

  constexpr bool usesOptimalTestBasis =
        is_OptimalTestSpace<
            typename std::tuple_element<std::tuple_size<TestSpaces>::value-1,
                                        TestSpaces>::type
        >::value;

  const unsigned int numS = x.size();

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  const size_t globalTotalTestSize
      = computeOffsets(globalTestSpaceOffsets, testSpaces);

  if(usesOptimalTestBasis)
  {
    for(size_t i=0; i<std::tuple_size<TestSpaces>::value; ++i)
    {
      globalTestSpaceOffsets[i] = 0;
    }
  }

  size_t globalTotalSolutionSize =
    computeOffsets(globalSolutionSpaceOffsets, solutionSpaces,
                   (!usesOptimalTestBasis)?globalTotalTestSize:0);
  globalTotalSolutionSize -= globalSolutionSpaceOffsets[0];

  const size_t globalSolutionSpaceOffset =
      globalSolutionSpaceOffsets[solutionSpaceIndex];

  scattering.resize(globalTotalTestSize
                    + ((!usesOptimalTestBasis)?globalTotalSolutionSize:0));
  scattering = 0;

  // Views on the FE bases on a single element
  auto testLocalViews     = getLocalViews(testSpaces);
  auto solutionLocalViews = getLocalViews(solutionSpaces);

  auto testLocalIndexSets     = getLocalIndexSets(testSpaces);
  auto solutionLocalIndexSets = getLocalIndexSets(solutionSpaces);

  for(const auto& e : elements(gridView)) {

    // Bind the local FE basis view to the current element
    /* TODO: only bind the space we use later */
    bindLocalViews(solutionLocalViews, e);
    bindLocalViews(testLocalViews, e);

    bindLocalIndexSets(solutionLocalIndexSets, solutionLocalViews);
    bindLocalIndexSets(testLocalIndexSets, testLocalViews);

    // Now get the local contribution to the right-hand side vector

    // Get the grid element from the local FE basis view
    typedef typename std::remove_reference<decltype(e)>::type Element;

    const int dim = Element::dimension;

    // Get set of shape functions for this element
    const auto& localFiniteElementTest =
        std::get<0>(testLocalViews).tree().finiteElement();
    const auto& localFiniteElementSolution =
        std::get<solutionSpaceIndex>(solutionLocalViews).tree().finiteElement();

    BlockVector<FieldVector<double,1> >
        localScattering(localFiniteElementTest.localBasis().size());
    localScattering = 0;

    /*TODO:
    - Find out what quadrature rule we are exactly using
    - Adapt quadrature to the kernel k also
    */
    int quadratureOrder = localFiniteElementSolution.localBasis().order()
                        + localFiniteElementTest.localBasis().order();
    /* Remark:
    - localFiniteElementTest.localBasis().order() is the degree of the polynomial
    - quad.size() is the number of quadrature points
    */

    const QuadratureRule<double, dim>& quad =
        QuadratureRules<double, dim>::rule(e.type(), quadratureOrder);

    // We get the values of the kernel at the quad points
    // and the angular directions. Result stored in kernelVect[iQuad][iS]
    std::vector<std::vector<double>> kernelVect;
    evaluateKernelElement(kernelVect, kernelS, e.geometry(), quad, sVector);

    // Loop over all quadrature points
    for ( size_t pt=0; pt < quad.size(); pt++ ) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();

      // The multiplicative factor in the integral transformation formula
      const double integrationElement = e.geometry().integrationElement(quadPos);

      // Evaluate all test shape function values at this quadrature point
      std::vector<FieldVector<double,1> > testShapeFunctionValues;
      localFiniteElementTest.localBasis()
          .evaluateFunction(quadPos, testShapeFunctionValues);

      for( unsigned int scatteringAngle=0;
           scatteringAngle<numS; ++scatteringAngle) {
        double uValue = 0; // in direction of scatteringAngle
        // Evaluate all shape function values at this point
        std::vector<FieldVector<double,1> > shapeFunctionValues;
        localFiniteElementSolution.localBasis().
            evaluateFunction(quadPos, shapeFunctionValues);
        for (size_t j=0; j<shapeFunctionValues.size(); j++)
        {
          /* This assumes that solutionLocalIndexSets and
           * globalSolutionSpaceOffset don't change for different
           * scattering angles.
           */
          auto row =
              std::get<solutionSpaceIndex>(solutionLocalIndexSets).index(j)[0]
            + globalSolutionSpaceOffset;
          uValue += x[scatteringAngle][row] * shapeFunctionValues[j];
        }

        const double factor = 1./numS * kernelVect[pt][scatteringAngle]
                              * uValue
                              * quad[pt].weight() * integrationElement;
        for (size_t i=0, i_max=localScattering.size(); i<i_max; i++)
          localScattering[i] += factor * testShapeFunctionValues[i];
      }

    }

    // TODO: We should probably not hard-code the test space index
    //       and testLocalOffset.
    auto rhsCopier
        = localToGlobalRHSCopier
            <typename std::remove_reference<decltype(localScattering)>::type,
             typename std::remove_reference<decltype(scattering)>::type,
             typename std::remove_reference<decltype(testLocalIndexSets)>::type,
             typename std::remove_reference<decltype(globalTestSpaceOffsets)>
                         ::type>
            (localScattering, scattering,
             testLocalIndexSets,
             {0},
             globalTestSpaceOffsets);
    rhsCopier(std::integral_constant<size_t, 0>{});

  }
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_SCATTERING_HH
