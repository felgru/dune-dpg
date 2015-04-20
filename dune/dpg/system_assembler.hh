// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_SYSTEM_ASSEMBLER_HH
#define DUNE_DPG_SYSTEM_ASSEMBLER_HH

#include <tuple>
#include <functional>
#include <memory>
#include <type_traits>

/* 7 would be enough, but better take some more, in case we
 * change our procedures later. */
#define BOOST_FUSION_INVOKE_PROCEDURE_MAX_ARITY 10

#include <boost/mpl/set.hpp>
#include <boost/mpl/transform.hpp>

#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/adapted/array.hpp>
#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/container/set/convert.hpp>
#include <boost/fusion/algorithm/auxiliary/copy.hpp>
#include <boost/fusion/algorithm/transformation/join.hpp>
#include <boost/fusion/algorithm/transformation/transform.hpp>
#include <boost/fusion/algorithm/transformation/zip.hpp>
#include <boost/fusion/algorithm/iteration/accumulate.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/functional/generation/make_fused_procedure.hpp>

#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dune/common/std/final.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include "assemble_helper.hh"
#include "assemble_types.hh"
#include "bilinearform.hh"
#include "innerproduct.hh"

namespace Dune {

// Compute the source term for a single element
template <class LocalViewTest, class LocalVolumeTerm>
void getVolumeTerm(const LocalViewTest& localViewTest,
                   BlockVector<FieldVector<double,1> >& localRhs,
                   LocalVolumeTerm&& localVolumeTerm)
{
  // Get the grid element from the local FE basis view
  typedef typename LocalViewTest::Element Element;
  const Element& element = localViewTest.element();

  const int dim = Element::dimension;

  // Get set of shape functions for this element
  const auto& localFiniteElementTest = localViewTest.tree().finiteElement();

  // Set all entries to zero
  localRhs.resize(localFiniteElementTest.localBasis().size());
  localRhs = 0;

  // A quadrature rule
  int order = dim*localFiniteElementTest.localBasis().order(); //TODO!!!!!!
  const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element.type(), order);


  // Loop over all quadrature points
  for ( size_t pt=0; pt < quad.size(); pt++ ) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = element.geometry().integrationElement(quadPos);

    double functionValue = localVolumeTerm(quadPos);

    // Evaluate all shape function values at this point
    std::vector<FieldVector<double,1> > shapeFunctionValues;
    localFiniteElementTest.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

    // Actually compute the vector entries
    for (size_t i=0; i<localRhs.size(); i++)
      localRhs[i] += shapeFunctionValues[i] * functionValue * quad[pt].weight() * integrationElement;

  }

}

struct getVolumeTermHelper
{
  template<class Seq>
  void operator()(const Seq& seq) const
  {
    using namespace boost::fusion;

    getVolumeTerm(*(at_c<0>(seq)), at_c<1>(seq), at_c<2>(seq));
    };
};

/**
 * class SystemAssembler
 *
 * \tparam TestSpaces     tuple of test spaces
 * \tparam SolutionSpaces tuple of solution spaces
 * \tparam BilinearForm   bilinear form describing the system
 * \tparam InnerProduct   inner product of the test space
 */
template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct>
class SystemAssembler
{
public:
  typedef typename boost::fusion::result_of::as_vector<
      typename boost::fusion::result_of::
      transform<TestSpaces, getLocalView>::type>::type TestLocalView;
  typedef typename boost::fusion::result_of::as_vector<
      typename boost::fusion::result_of::
      transform<SolutionSpaces, getLocalView>::type>::type SolutionLocalView;

  SystemAssembler () = delete;
  constexpr SystemAssembler (TestSpaces     testSpaces,
                             SolutionSpaces solutionSpaces,
                             BilinearForm   bilinearForm,
                             InnerProduct   innerProduct)
             : testSpaces(testSpaces),
               solutionSpaces(solutionSpaces),
               bilinearForm(bilinearForm),
               innerProduct(innerProduct)
  { };

  template <class VolumeTerms>
  void assembleSystem(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                      BlockVector<FieldVector<double,1> >& rhs,
                      VolumeTerms&& volumeTerms);

  template <SpaceType spaceType, size_t spaceIndex, class ValueType>
  void applyDirichletBoundary(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                              BlockVector<FieldVector<double,1> >& rhs,
                              const std::vector<bool>& dirichletNodes,
                              const ValueType& value);

  const TestSpaces& getTestSpaces() const
  { return testSpaces; };

  const SolutionSpaces& getSolutionSpaces() const
  { return solutionSpaces; };

private:
  TestSpaces     testSpaces;
  SolutionSpaces solutionSpaces;
  BilinearForm   bilinearForm;
  InnerProduct   innerProduct;
};

template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct>
auto make_SystemAssembler(TestSpaces     testSpaces,
                          SolutionSpaces solutionSpaces,
                          BilinearForm   bilinearForm,
                          InnerProduct   innerProduct)
    -> SystemAssembler<TestSpaces, SolutionSpaces,
                       BilinearForm, InnerProduct>
{
  return SystemAssembler<TestSpaces, SolutionSpaces,
                         BilinearForm, InnerProduct>
                      (testSpaces,
                       solutionSpaces,
                       bilinearForm,
                       innerProduct);
}

template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct>
template <class VolumeTerms>
void SystemAssembler<TestSpaces, SolutionSpaces, BilinearForm, InnerProduct>::
assembleSystem(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
               BlockVector<FieldVector<double,1> >& rhs,
               VolumeTerms&& volumeTerms)
{
  using namespace boost::fusion;

  // Get the grid view from the finite element basis
  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  /* TODO: make this a vector of pointers to localVolumeTerm */
  auto localVolumeTerms =
      as_vector(transform(volumeTerms,
                          getLocalVolumeTerm<GridView>(gridView)));

  auto testBasisIndexSet = as_vector(transform(testSpaces, getIndexSet()));
  auto solutionBasisIndexSet = as_vector(transform(solutionSpaces,
                                                   getIndexSet()));

  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  fold(zip(globalTestSpaceOffsets, testBasisIndexSet), 0, globalOffsetHelper());
  size_t globalTotalTestSize =
      globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value-1]
      + at_c<std::tuple_size<TestSpaces>::value-1>(testBasisIndexSet).size();

  fold(zip(globalSolutionSpaceOffsets, solutionBasisIndexSet),
       globalTotalTestSize, globalOffsetHelper());
  size_t globalTotalSolutionSize =
      globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value-1]
      + at_c<std::tuple_size<SolutionSpaces>::value-1>
            (solutionBasisIndexSet).size()
      - globalTotalTestSize;

  auto n = globalTotalTestSize + globalTotalSolutionSize;

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  MatrixIndexSet occupationPattern;
  occupationPattern.resize(n, n);
  bilinearForm.getOccupationPattern(occupationPattern);
  innerProduct.getOccupationPattern(occupationPattern);

  /* Add the diagonal of the matrix, since we need it for the
   * interpolation of boundary values. */
  for (size_t i=0; i<n; i++)
  {
    occupationPattern.add(i, i);
  }

  // ... and give it the occupation pattern we want.
  occupationPattern.exportIdx(matrix);

  // set rhs to correct length -- the total number of basis vectors in the basis
  rhs.resize(globalTotalTestSize + globalTotalSolutionSize);

  // Set all entries to zero
  matrix = 0;
  rhs = 0;

  // Views on the FE bases on a single element
  auto solutionLocalView = as_vector(transform(solutionSpaces, getLocalView()));
  auto testLocalView     = as_vector(transform(testSpaces, getLocalView()));

  auto solutionLocalIndexSet = as_vector(transform(solutionBasisIndexSet,
                                                   getLocalIndexSet()));
  auto testLocalIndexSet     = as_vector(transform(testBasisIndexSet,
                                                   getLocalIndexSet()));

  /* create sets of index pairs to loop over. */
  typedef typename boost::mpl::fold<
      typename boost::mpl::transform<
          /* This as_vector is probably not needed for boost::fusion 1.58
           * or higher. */
          typename result_of::as_vector<typename std::remove_reference<
                decltype(bilinearForm.getTerms())>::type
              >::type
        , firstTwo::mpl<boost::mpl::_1>
        >::type
    , boost::mpl::set0<>
    , boost::mpl::insert<boost::mpl::_1,boost::mpl::_2>
    >::type BFIndices;
  typedef typename boost::mpl::fold<
      typename boost::mpl::transform<
          typename result_of::as_vector<typename std::remove_reference<
                decltype(innerProduct.getTerms())>::type
              >::type
        , firstTwo::mpl<boost::mpl::_1>
        >::type
    , boost::mpl::set0<>
    , boost::mpl::insert<boost::mpl::_1,boost::mpl::_2>
    >::type IPIndices;

  auto bfIndices = BFIndices{};
  auto ipIndices = IPIndices{};

  // A loop over all elements of the grid
  for(const auto& e : elements(gridView)) {

    // Bind the local FE basis view to the current element
    for_each(solutionLocalView, applyBind<decltype(e)>(e));
    for_each(testLocalView, applyBind<decltype(e)>(e));

    for_each(zip(solutionLocalIndexSet, solutionLocalView),
             make_fused_procedure(bindLocalIndexSet()));
    for_each(zip(testLocalIndexSet, testLocalView),
             make_fused_procedure(bindLocalIndexSet()));

    /* Bind the bilinearForm and the innerProduct to the local views. */
    bilinearForm.bind(testLocalView, solutionLocalView);
    innerProduct.bind(testLocalView);

    // Now let's get the element stiffness matrix and the Gram matrix
    // for the test space.
    Matrix<FieldMatrix<double,1,1> > bfElementMatrix;
    Matrix<FieldMatrix<double,1,1> > ipElementMatrix;

    bilinearForm.getLocalMatrix(bfElementMatrix);
    innerProduct.getLocalMatrix(ipElementMatrix);


    // Add element stiffness matrix onto the global stiffness matrix
    auto cp = fused_procedure<localToGlobalCopier<decltype(ipElementMatrix),
                        typename remove_reference<decltype(matrix)>::type> >
                    (localToGlobalCopier<decltype(ipElementMatrix),
                        typename remove_reference<decltype(matrix)>::type>
                                        (ipElementMatrix, matrix));
    auto cpm = fused_procedure<localToGlobalCopier<decltype(bfElementMatrix),
                        typename remove_reference<decltype(matrix)>::type,
                        true> >
                    (localToGlobalCopier<decltype(bfElementMatrix),
                        typename remove_reference<decltype(matrix)>::type,
                        true>
                                        (bfElementMatrix, matrix));

    /* iterate over bfIndices and ipIndices */
    auto testZip = zip(testLocalView,
                       testLocalIndexSet,
                       bilinearForm.getLocalTestSpaceOffsets(),
                       globalTestSpaceOffsets);
    auto solutionZip = zip(solutionLocalView,
                           solutionLocalIndexSet,
                           bilinearForm.getLocalSolutionSpaceOffsets(),
                           globalSolutionSpaceOffsets);
    for_each(ipIndices,
             localToGlobalCopyHelper<decltype(testZip),
                                     decltype(testZip),
                                     decltype(cp)>
                                    (testZip, testZip, std::ref(cp)));
    for_each(bfIndices,
             localToGlobalCopyHelper<decltype(solutionZip),
                                     decltype(testZip),
                                     decltype(cpm)>
                                    (solutionZip, testZip, std::ref(cpm)));

    /* TODO: Make this work for larger sets of test spaces. */
    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> >
        localRhs[std::tuple_size<
                 typename std::remove_reference<VolumeTerms>::type>::value];
    for_each(localVolumeTerms, applyBind<decltype(e)>(e));

    using RHSZipHelper = vector<decltype(testLocalView)&,
                                decltype(localRhs)&,
                                decltype(localVolumeTerms)&>;
    for_each(zip_view<RHSZipHelper>(RHSZipHelper(testLocalView,
                                                 localRhs,
                                                 localVolumeTerms)),
             getVolumeTermHelper());

    auto cpr = fused_procedure<localToGlobalRHSCopier<
                        typename remove_reference<decltype(rhs)>::type> >
                    (localToGlobalRHSCopier<
                        typename remove_reference<decltype(rhs)>::type>(rhs));
    for_each(zip(localRhs,
                 testLocalIndexSet,
                 globalTestSpaceOffsets),
             std::ref(cpr));

  }

  /* free memory handled by raw pointers */
  for_each(testLocalIndexSet,     default_deleter());
  for_each(solutionLocalIndexSet, default_deleter());
  for_each(testLocalView,         default_deleter());
  for_each(solutionLocalView,     default_deleter());
}

template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct>
template <SpaceType spaceType, size_t spaceIndex, class ValueType>
void SystemAssembler<TestSpaces, SolutionSpaces, BilinearForm, InnerProduct>::
applyDirichletBoundary(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                       BlockVector<FieldVector<double,1> >& rhs,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& value)
{
  using namespace boost::fusion;
  static_assert(std::is_arithmetic<ValueType>::value,
                "applyDirichletBoundary not implemented for non arithmetic "
                "boundary data types.");

  size_t spaceSize;
  if(spaceType==SpaceType::test)
    spaceSize = std::get<spaceIndex>(testSpaces).indexSet().size();
  else
    spaceSize = std::get<spaceIndex>(solutionSpaces).indexSet().size();

  size_t globalOffset;
  {
    // Total number of degrees of freedom
    auto testBasisIndexSet = as_vector(transform(testSpaces, getIndexSet()));
    auto solutionBasisIndexSet = as_vector(transform(solutionSpaces,
                                                     getIndexSet()));

    /* set up global offsets */
    size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
    size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
    fold(zip(globalTestSpaceOffsets, testBasisIndexSet),
         0, globalOffsetHelper());
    size_t globalTotalTestSize =
        globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value-1]
        + at_c<std::tuple_size<TestSpaces>::value-1>(testBasisIndexSet).size();

    fold(zip(globalSolutionSpaceOffsets, solutionBasisIndexSet),
         globalTotalTestSize, globalOffsetHelper());

    if(spaceType==SpaceType::test)
      globalOffset = globalTestSpaceOffsets[spaceIndex];
    else
      globalOffset = globalSolutionSpaceOffsets[spaceIndex];
  }

  // Set Dirichlet values
  for (size_t i=0; i<spaceSize; i++)
  {
    if (dirichletNodes[i])
    {
      /* TODO: Needs adaptation when value is a function. */
      rhs[globalOffset+i] = value;
    }
  }

  ////////////////////////////////////////////
  //   Modify Dirichlet rows
  ////////////////////////////////////////////

  // loop over the matrix rows
  for (size_t i=0; i<spaceSize; i++)
  {
    if (dirichletNodes[i])
    {
      auto cIt    = matrix[globalOffset+i].begin();
      auto cEndIt = matrix[globalOffset+i].end();
      // loop over nonzero matrix entries in current row
      for (; cIt!=cEndIt; ++cIt)
      {
        if (globalOffset+i==cIt.index())
        {
          *cIt = 1;
        }
        else
        {
          /* TODO: This seems somewhat redundant... */
          *cIt = 0;
          matrix[cIt.index()][globalOffset+i]=0;
        }
      }
    }

  }
};

} // end namespace Dune

#endif // DUNE_DPG_SYSTEM_ASSEMBLER_HH
