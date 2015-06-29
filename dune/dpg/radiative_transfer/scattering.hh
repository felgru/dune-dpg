// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_SCATTERING_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_SCATTERING_HH

#include <tuple>
#include <functional>
#include <memory>
#include <type_traits>

#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/assemble_helper.hh>

#include <boost/fusion/algorithm/iteration/for_each.hpp>

#include <boost/math/constants/constants.hpp>

namespace Dune {

/**
 * \brief This constructs the right hand side vector of a DPG system.
 *
 * \tparam TestSpaces      tuple of test spaces
 * \tparam SolutionSpaces  tuple of solution spaces
 * \tparam FormulationType either SaddlepointFormulation or DPGFormulation
 */
template<class TestSpaces,
         class SolutionSpaces,
         class FormulationType>
class ScatteringAssembler
{
public:
  ScatteringAssembler () = delete;
  /**
   * \brief constructor for ScatteringAssembler
   *
   * \note For your convenience, use make_ScatteringAssembler() instead.
   */
  constexpr ScatteringAssembler (TestSpaces testSpaces,
                                 SolutionSpaces solutionSpaces)
             : testSpaces(testSpaces),
               solutionSpaces(solutionSpaces)
  { };

  /**
   * \brief Assemble the vector corresponding to the scattering integral
   * for a given set of discrete ordinate solutions.
   *
   * \todo: For now, the scattering kernel is assumed to be constant 1.
   *
   * \param[out] rhs  the scattering vector
   * \param[in]  i    index of the scattering direction
   * \param[in]  x    the vectors of the solutions of the previous iteration
   */
  template<size_t solutionSpaceIndex>
  void assembleScattering
         (BlockVector<FieldVector<double,1> >& rhs,
          int i,
          const std::vector<BlockVector<FieldVector<double,1> >>& x);

private:
  TestSpaces     testSpaces;
  SolutionSpaces solutionSpaces;
};

/**
 * \brief Creates an ScatteringAssembler,
 *        deducing the target type from the types of arguments.
 *
 * \param testSpaces     a tuple of test spaces
 * \param testSpaces     a tuple of solution spaces
 * \tparam FormulationType either SaddlepointFormulation or DPGFormulation
 */
template<class TestSpaces,
         class SolutionSpaces,
         class FormulationType>
auto make_ScatteringAssembler(TestSpaces testSpaces,
                              SolutionSpaces solutionSpaces,
                              FormulationType)
     -> ScatteringAssembler<TestSpaces, SolutionSpaces, FormulationType>
{
  return ScatteringAssembler<TestSpaces, SolutionSpaces,
                             FormulationType>
                            (testSpaces, solutionSpaces);
}

template<class TestSpaces,
         class SolutionSpaces,
         class FormulationType>
template<size_t solutionSpaceIndex>
void ScatteringAssembler<TestSpaces, SolutionSpaces, FormulationType>::
assembleScattering(BlockVector<FieldVector<double,1> >& rhs,
                   int i,
                   const std::vector<BlockVector<FieldVector<double,1> >>& x)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  constexpr bool isSaddlepoint =
        std::is_same<
             typename std::decay<FormulationType>::type
           , SaddlepointFormulation
        >::value;

  const int numS = x.size();

  // Get the grid view from the finite element basis
  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  auto testBasisIndexSet     = as_vector(transform(testSpaces,
                                                   getIndexSet()));
  auto solutionBasisIndexSet = as_vector(transform(solutionSpaces,
                                                   getIndexSet()));

  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  size_t globalTotalTestSize = 0;

  if(isSaddlepoint) {
    globalTotalTestSize =
      fold(zip(globalTestSpaceOffsets, testBasisIndexSet),
           0, globalOffsetHelper());
  } else { /* DPG formulation */
    for(size_t i=0; i<std::tuple_size<TestSpaces>::value; ++i)
    {
      globalTestSpaceOffsets[i] = 0;
    }
  }

  size_t globalTotalSolutionSize =
    fold(zip(globalSolutionSpaceOffsets, solutionBasisIndexSet),
         isSaddlepoint?globalTotalTestSize:0, globalOffsetHelper());
  globalTotalSolutionSize -= globalSolutionSpaceOffsets[0];

  double globalSolutionSpaceOffset =
      globalSolutionSpaceOffsets[solutionSpaceIndex];

  if(!isSaddlepoint) globalTotalTestSize = globalTotalSolutionSize;

  // set rhs to the right size
  rhs.resize(globalTotalTestSize
             + (isSaddlepoint?globalTotalSolutionSize:0));
  rhs = 0;

  // Views on the FE bases on a single element
  auto testLocalView     = as_vector(transform(testSpaces, getLocalView()));
  auto solutionLocalView = as_vector(transform(solutionSpaces, getLocalView()));

  auto testLocalIndexSet     = as_vector(transform(testBasisIndexSet,
                                                   getLocalIndexSet()));
  auto solutionLocalIndexSet = as_vector(transform(solutionBasisIndexSet,
                                                   getLocalIndexSet()));

  // A loop over all elements of the grid
  for(const auto& e : elements(gridView)) {

    // Bind the local FE basis view to the current element
    /* TODO: only bind the space we use later */
    for_each(solutionLocalView, applyBind<decltype(e)>(e));
    for_each(testLocalView, applyBind<decltype(e)>(e));

    for_each(zip(solutionLocalIndexSet, solutionLocalView),
             make_fused_procedure(bindLocalIndexSet()));
    for_each(zip(testLocalIndexSet, testLocalView),
             make_fused_procedure(bindLocalIndexSet()));

    // Now get the local contribution to the right-hand side vector

    // Get the grid element from the local FE basis view
    typedef typename std::remove_reference<decltype(e)>::type Element;

    const int dim = Element::dimension;

    // Get set of shape functions for this element
    const auto& localFiniteElementTest =
        at_c<0>(testLocalView)->tree().finiteElement();
    const auto& localFiniteElementSolution =
        at_c<solutionSpaceIndex>(solutionLocalView)->tree().finiteElement();

    BlockVector<FieldVector<double,1> >
        localScattering(localFiniteElementTest.localBasis().size());

    // Set all entries to zero
    localScattering = 0;

    // A quadrature rule
    int order = dim*localFiniteElementTest.localBasis().order(); //TODO!!!!!!
    const QuadratureRule<double, dim>& quad =
        QuadratureRules<double, dim>::rule(e.type(), order);


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

      for( size_t scatteringAngle=0;
           scatteringAngle<numS; ++scatteringAngle) {
        double uValue = 0; // in direction of scatteringAngle
        // Evaluate all shape function values at this point
        std::vector<FieldVector<double,1> > shapeFunctionValues;
        localFiniteElementSolution.localBasis().
            evaluateFunction(quadPos, shapeFunctionValues);
        for (size_t j=0; j<shapeFunctionValues.size(); j++)
        {
          /* This assumes that solutionLocalIndexSet and
           * globalSolutionSpaceOffset don't change for different
           * scattering angles.
           */
          auto row =
              at_c<solutionSpaceIndex>(solutionLocalIndexSet)->index(j)[0]
            + globalSolutionSpaceOffset;
          uValue += x[scatteringAngle][row] * shapeFunctionValues[j];
        }

        // Actually compute the vector entries
        const double factor = 1./numS
                              * uValue * quad[pt].weight() * integrationElement;
        for (size_t i=0, i_max=localScattering.size(); i<i_max; i++)
          localScattering[i] += factor * testShapeFunctionValues[i];
      }

    }

    auto rhsCopier = localToGlobalRHSCopier
                     <typename remove_reference<decltype(rhs)>::type>(rhs);
    rhsCopier(localScattering,
              at_c<0>(testLocalIndexSet),
              globalTestSpaceOffsets[0]
             );

  }

  /* free memory handled by raw pointers */
  for_each(solutionLocalIndexSet, default_deleter());
  for_each(testLocalIndexSet,     default_deleter());
  for_each(solutionLocalView,     default_deleter());
  for_each(testLocalView,         default_deleter());
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_SCATTERING_HH
