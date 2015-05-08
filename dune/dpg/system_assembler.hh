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

#include <boost/mpl/identity.hpp>
#include <boost/mpl/range_c.hpp>
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

namespace detail {

template<class Form, class newTestSpace, class newFormulationType>
struct replaceTestSpaceAndType {};

template<class TestSpaces, class SolutionSpaces, class BilinearTerms,
         class FormulationType, class newTestSpaces, class newFormulationType>
struct replaceTestSpaceAndType<BilinearForm<TestSpaces, SolutionSpaces,
                                            BilinearTerms, FormulationType>,
                               newTestSpaces, newFormulationType>
{
    typedef BilinearForm<newTestSpaces, SolutionSpaces,
                         BilinearTerms, newFormulationType> type;
};

template<class Form, class newTestSpace>
struct replaceTestSpace {};

template<class TestSpaces, class InnerProductTerms, class newTestSpaces>
struct replaceTestSpace<InnerProduct<TestSpaces, InnerProductTerms>,
                        newTestSpaces>
{
    typedef InnerProduct<newTestSpaces, InnerProductTerms> type;
};

} // end namespace detail

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

namespace detail {
struct getVolumeTermHelper
{
  template<class Seq>
  void operator()(const Seq& seq) const
  {
    using namespace boost::fusion;

    getVolumeTerm(*(at_c<0>(seq)), at_c<1>(seq), at_c<2>(seq));
    };
};
} // end namespace detail

/**
 * \brief This constructs the matrix and vector of a DPG system.
 *
 * \tparam TestSpaces     tuple of test spaces
 * \tparam SolutionSpaces tuple of solution spaces
 * \tparam BilinForm      bilinear form describing the system
 * \tparam InProduct      inner product of the test space
 * \tparam FormulationType either SaddlepointFormulation or DPGFormulation
 */
template<class TestSpaces, class SolutionSpaces,
         class BilinForm, class InProduct,
         class FormulationType>
class SystemAssembler
{
public:
  //! tuple type for the local views of the test spaces
  typedef typename boost::fusion::result_of::as_vector<
      typename boost::fusion::result_of::
      transform<TestSpaces, detail::getLocalView>::type>::type TestLocalView;
  //! tuple type for the local views of the solution spaces
  typedef typename boost::fusion::result_of::as_vector<
      typename boost::fusion::result_of::
      transform<SolutionSpaces, detail::getLocalView>::type
      >::type SolutionLocalView;
  //! type of the bilinear form describing this DPG system
  typedef typename std::conditional<
        std::is_same<
             typename std::decay<FormulationType>::type
           , SaddlepointFormulation
        >::value
      , BilinForm
      , typename detail::replaceTestSpaceAndType<
                          BilinForm, TestSpaces, FormulationType
                         >::type
      >::type BilinearForm;
  //! type of the inner product on the test spaces
  typedef typename std::conditional<
        std::is_same<
             typename std::decay<FormulationType>::type
           , SaddlepointFormulation
        >::value
      , InProduct
      , typename detail::replaceTestSpace<InProduct, TestSpaces>::type
      >::type InnerProduct;

  SystemAssembler () = delete;
  /**
   * \brief constructor for SystemAssembler
   *
   * \note For your convenience, use make_SystemAssembler() instead.
   */
  constexpr SystemAssembler (TestSpaces     testSpaces,
                             SolutionSpaces solutionSpaces,
                             BilinForm      bilinearForm,
                             InProduct      innerProduct)
             : testSpaces(testSpaces),
               solutionSpaces(solutionSpaces),
               bilinearForm(detail::make_BilinearForm<FormulationType>(testSpaces,
                                              solutionSpaces,
                                              bilinearForm.getTerms())),
               innerProduct(make_InnerProduct(testSpaces,
                                              innerProduct.getTerms()))
  { };

  /**
   * \brief Assemble the DPG system for a given rhs function.
   *
   * Given a tuple of right hand side functions \p volumeTerms,
   * this assembles the matrix and vector of the corresponding
   * DPG system. \p matrix and \p rhs will be overwritten by this
   * function.
   *
   * \param[out] matrix      the matrix of the DPG system
   * \param[out] rhs         the rhs vector of the DPG system
   * \param[in]  volumeTerms the rhs functions describing the DPG system
   * \tparam     VolumeTerms a tuple type of rhs functions
   */
  template <class VolumeTerms>
  void assembleSystem(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                      BlockVector<FieldVector<double,1> >& rhs,
                      VolumeTerms&& volumeTerms);

  /**
   * \brief Apply Dirichlet boundary values on a test space
   *
   * \param[in,out] matrix      the matrix of the DPG system
   * \param[in,out] rhs         the rhs vector of the DPG system
   * \param[in] dirichletNodes  true marks the dofs in the Dirichlet boundary
   * \param[in] value           the Dirichlet boundary value
   * \tparam spaceIndex  the index of the test space on which we apply
   *                     the boundary data
   * \tparam ValueType   we take either constants or functions for \p value
   */
  template <size_t spaceIndex, class ValueType>
  void applyDirichletBoundaryTest(
                              BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                              BlockVector<FieldVector<double,1> >& rhs,
                              const std::vector<bool>& dirichletNodes,
                              const ValueType& value);

  /**
   * \brief Apply Dirichlet boundary values on a solution space
   *
   * \param[in,out] matrix      the matrix of the DPG system
   * \param[in,out] rhs         the rhs vector of the DPG system
   * \param[in] dirichletNodes  true marks the dofs in the Dirichlet boundary
   * \param[in] value           the Dirichlet boundary value
   * \tparam spaceIndex  the index of the solution space on which we apply
   *                     the boundary data
   * \tparam ValueType   we take either constants or functions for \p value
   */
  template <size_t spaceIndex, class ValueType>
  void applyDirichletBoundarySolution(
                              BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                              BlockVector<FieldVector<double,1> >& rhs,
                              const std::vector<bool>& dirichletNodes,
                              const ValueType& value);

  template <size_t spaceIndex, unsigned int dim>
  void applyWeakBoundaryCondition(
                              BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                              FieldVector<double, dim> beta,
                              double c);

  /**
   * \brief Does exactly what it says on the tin.
   */
  const TestSpaces& getTestSpaces() const
  { return testSpaces; };

  /**
   * \brief Does exactly what it says on the tin.
   */
  const SolutionSpaces& getSolutionSpaces() const
  { return solutionSpaces; };

private:

  template <SpaceType spaceType, size_t spaceIndex,
            class ValueType, class Spaces>
  void applyDirichletBoundaryImpl(
                       BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                       BlockVector<FieldVector<double,1> >& rhs,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& value,
                       const Spaces& spaces);

  TestSpaces     testSpaces;
  SolutionSpaces solutionSpaces;
  BilinearForm   bilinearForm;
  InnerProduct   innerProduct;
};

/**
 * \brief Creates a SystemAssembler,
 *        deducing the target type from the types of arguments.
 *
 * \param testSpaces     a tuple of test spaces
 * \param solutionSpaces a tuple of solution spaces
 * \param bilinearForm   the bilinear form describing the DPG system
 * \param innerProduct   the inner product of the test spaces
 * \tparam FormulationType either SaddlepointFormulation or DPGFormulation
 */
template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct,
         class FormulationType>
auto make_SystemAssembler(TestSpaces     testSpaces,
                          SolutionSpaces solutionSpaces,
                          BilinearForm   bilinearForm,
                          InnerProduct   innerProduct,
                          FormulationType)
    -> SystemAssembler<TestSpaces, SolutionSpaces,
                       BilinearForm, InnerProduct,
                       FormulationType>
{
  return SystemAssembler<TestSpaces, SolutionSpaces,
                         BilinearForm, InnerProduct,
                         FormulationType>
                      (testSpaces,
                       solutionSpaces,
                       bilinearForm,
                       innerProduct);
}

template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct,
         class FormulationType>
template <class VolumeTerms>
void SystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct, FormulationType>::
assembleSystem(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
               BlockVector<FieldVector<double,1> >& rhs,
               VolumeTerms&& volumeTerms)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  constexpr bool isSaddlepoint =
        std::is_same<
             typename std::decay<FormulationType>::type
           , SaddlepointFormulation
        >::value;

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
  size_t globalTotalTestSize = 0;

  if(isSaddlepoint) {
    fold(zip(globalTestSpaceOffsets, testBasisIndexSet),
         0, globalOffsetHelper());
    globalTotalTestSize =
        globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value-1]
        + at_c<std::tuple_size<TestSpaces>::value-1>(testBasisIndexSet).size();
  } else { /* DPG formulation */
    for(size_t i=0; i<std::tuple_size<TestSpaces>::value; ++i)
    {
      globalTestSpaceOffsets[i] = 0;
    }
  }

  fold(zip(globalSolutionSpaceOffsets, solutionBasisIndexSet),
       isSaddlepoint?globalTotalTestSize:0, globalOffsetHelper());

  size_t globalTotalSolutionSize =
      globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value-1]
      + at_c<std::tuple_size<SolutionSpaces>::value-1>
            (solutionBasisIndexSet).size()
      - globalSolutionSpaceOffsets[0];

  auto n = globalTotalSolutionSize;
  if(isSaddlepoint)
  {
    n+=globalTotalTestSize;
  }

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  MatrixIndexSet occupationPattern;
  occupationPattern.resize(n, n);
  bilinearForm.template getOccupationPattern<isSaddlepoint>
               (occupationPattern,
                0, isSaddlepoint?globalTotalTestSize:0);
  if(isSaddlepoint)
  {
    innerProduct.getOccupationPattern(occupationPattern);

    /* Add the diagonal of the matrix, since we need it for the
     * interpolation of boundary values. */
    for (size_t i=0; i<n; i++)
    {
      occupationPattern.add(i, i);
    }
  }

  // ... and give it the occupation pattern we want.
  occupationPattern.exportIdx(matrix);

  // set rhs to correct length -- the total number of basis vectors in the basis
  rhs.resize((isSaddlepoint?globalTotalTestSize:0) + globalTotalSolutionSize);

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
        , mpl::firstTwo<boost::mpl::_1>
        >::type
    , boost::mpl::set0<>
    , boost::mpl::insert<boost::mpl::_1,boost::mpl::_2>
    >::type BFIndices;
  typedef typename boost::mpl::fold<
      typename boost::mpl::transform<
          typename result_of::as_vector<typename std::remove_reference<
                decltype(innerProduct.getTerms())>::type
              >::type
        , mpl::firstTwo<boost::mpl::_1>
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
    if(isSaddlepoint)
    {
      innerProduct.bind(testLocalView);
    }

    // Now let's get the element stiffness matrix and the Gram matrix
    // for the test space.
    Matrix<FieldMatrix<double,1,1> > bfElementMatrix;
    Matrix<FieldMatrix<double,1,1> > ipElementMatrix;

    bilinearForm.getLocalMatrix(bfElementMatrix);
    if(isSaddlepoint) {
      innerProduct.getLocalMatrix(ipElementMatrix);
    }


    // Add element stiffness matrix onto the global stiffness matrix
    auto cp = fused_procedure<localToGlobalCopier<decltype(ipElementMatrix),
                        typename remove_reference<decltype(matrix)>::type> >
                    (localToGlobalCopier<decltype(ipElementMatrix),
                        typename remove_reference<decltype(matrix)>::type>
                                        (ipElementMatrix, matrix));
    auto cpm = fused_procedure<localToGlobalCopier<decltype(bfElementMatrix),
                        typename remove_reference<decltype(matrix)>::type,
                        isSaddlepoint> >
                    (localToGlobalCopier<decltype(bfElementMatrix),
                        typename remove_reference<decltype(matrix)>::type,
                        isSaddlepoint>
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
    if(isSaddlepoint)
    {
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
    } else {
      typedef
          typename boost::mpl::transform<
              typename result_of::as_vector<typename boost::mpl::range_c<
                                  size_t,0,std::tuple_size<SolutionSpaces
                               >::value>::type
                  >::type
            , mpl::tupleOf0And<boost::mpl::_1>
            >::type BFIndices;

      auto bfIndices = BFIndices{};
      for_each(bfIndices,
               localToGlobalCopyHelper<decltype(solutionZip),
                                       decltype(testZip),
                                       decltype(cpm)>
                                      (solutionZip, testZip, std::ref(cpm)));
    }

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

    /* TODO: This will break with more than 1 test space having a rhs! */
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
         class BilinearForm, class InnerProduct,
         class FormulationType>
template <size_t spaceIndex, class ValueType>
void SystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct, FormulationType>::
applyDirichletBoundaryTest
                      (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                       BlockVector<FieldVector<double,1> >& rhs,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& value)
{
  applyDirichletBoundaryImpl<SpaceType::test, spaceIndex,
                             ValueType, TestSpaces>
                             (matrix, rhs, dirichletNodes,
                              value, testSpaces);
}

template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct,
         class FormulationType>
template <size_t spaceIndex, class ValueType>
void SystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct, FormulationType>::
applyDirichletBoundarySolution
                      (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                       BlockVector<FieldVector<double,1> >& rhs,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& value)
{
  applyDirichletBoundaryImpl<SpaceType::solution, spaceIndex,
                             ValueType, SolutionSpaces>
                             (matrix, rhs, dirichletNodes,
                              value, solutionSpaces);
}

template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct,
         class FormulationType>
template <SpaceType spaceType, size_t spaceIndex, class ValueType, class Spaces>
void SystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct, FormulationType>::
applyDirichletBoundaryImpl(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                       BlockVector<FieldVector<double,1> >& rhs,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& value,
                       const Spaces& spaces)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  static_assert(std::is_arithmetic<ValueType>::value,
                "applyDirichletBoundary not implemented for non arithmetic "
                "boundary data types.");

  constexpr bool isSaddlepoint =
        std::is_same<
             typename std::decay<FormulationType>::type
           , SaddlepointFormulation
        >::value;

  const size_t spaceSize =
        std::get<spaceIndex>(spaces).indexSet().size();

  size_t globalOffset;
  {
    // Total number of degrees of freedom
    auto testBasisIndexSet = as_vector(transform(testSpaces, getIndexSet()));
    auto solutionBasisIndexSet = as_vector(transform(solutionSpaces,
                                                     getIndexSet()));

    /* set up global offsets */
    size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
    size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
    size_t globalTotalTestSize = 0;

    if(isSaddlepoint) {
      fold(zip(globalTestSpaceOffsets, testBasisIndexSet),
           0, globalOffsetHelper());
      globalTotalTestSize =
          globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value-1]
          + at_c<std::tuple_size<TestSpaces>::value-1>(testBasisIndexSet).size();
    } else { /* DPG formulation */
      for(size_t i=0; i<std::tuple_size<TestSpaces>::value; ++i)
      {
        globalTestSpaceOffsets[i] = 0;
      }
    }

    fold(zip(globalSolutionSpaceOffsets, solutionBasisIndexSet),
         isSaddlepoint?globalTotalTestSize:0, globalOffsetHelper());

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

template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct,
         class FormulationType>
template <size_t spaceIndex, unsigned int dim>
void SystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct, FormulationType>::
applyWeakBoundaryCondition
                    (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                     FieldVector<double, dim> beta,
                     double c)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  static_assert(std::is_same<
             typename std::decay<FormulationType>::type
           , DPGFormulation
                            >::value,
                "applyWeakBoundaryConditions not implemented "
                "for Saddlepointformulation ");

  // Get the grid view from the finite element basis
  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  // get global solution Space Offsets
  auto solutionBasisIndexSet = as_vector(transform(solutionSpaces,
                                                     getIndexSet()));

  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  fold(zip(globalSolutionSpaceOffsets, solutionBasisIndexSet),
       0, globalOffsetHelper());
  size_t globalOffset = globalSolutionSpaceOffsets[spaceIndex];
  // get local view and local index set for the spaceIndex'th solution space
  auto localView = std::get<spaceIndex>(solutionSpaces).localView();
  typedef decltype(localView) LocalView;
  auto localIndexSet = at_c<spaceIndex>(solutionBasisIndexSet).localIndexSet();

  for(const auto& e : elements(gridView))
  {
    localView.bind(e);
    localIndexSet.bind(localView);

    const auto& localFiniteElement = localView.tree().finiteElement();

    int order = 2*(dim*localFiniteElement.localBasis().order()-1);

    size_t n = localFiniteElement.localBasis().size();

    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    // Set all matrix entries to zero
    elementMatrix.setSize(n,n);
    elementMatrix = 0;      // fills the entire matrix with zeroes

    for (auto&& intersection : intersections(gridView, e))
    {
      if (intersection.boundary()) //if the intersection is at the (physical) boundary of the domain
      {
        const FieldVector<double,dim>& centerOuterNormal =
                intersection.centerUnitOuterNormal();

        if ((beta*centerOuterNormal) > -1e-10) //everywhere except inflow boundary
        {
          const QuadratureRule<double, dim-1>& quadFace =
                  QuadratureRules<double, dim-1>::rule(intersection.type(), order);

          for (size_t pt=0; pt < quadFace.size(); pt++)
          {
            // Position of the current quadrature point in the reference element (face!)
            const FieldVector<double,dim-1>& quadFacePos = quadFace[pt].position();

            const double integrationElement = intersection.geometry().integrationElement(quadFacePos);

            // position of the quadrature point within the element
            const FieldVector<double,dim> elementQuadPos =
                intersection.geometryInInside().global(quadFacePos);

            // values of the shape functions
            std::vector<FieldVector<double,1> > solutionValues;
            localFiniteElement.localBasis().evaluateFunction(elementQuadPos,
                                                        solutionValues);
            for (size_t i=0; i<n; i++)
            {
              for (size_t j=0; j<n; j++)
              {
                elementMatrix[i][j]
                        += (c * solutionValues[i] * solutionValues[j])
                           * quadFace[pt].weight() * integrationElement;
              }
            }

          }
        }
      }
    }
    for (size_t i=0; i<n; i++)
    {
      auto row = localIndexSet.index(i)[0];
      for (size_t j=0; j<n; j++)
      {
        auto col = localIndexSet.index(j)[0];
        matrix[row+globalOffset][col+globalOffset]
                        += elementMatrix[i][j];

      }
    }
  }
};


} // end namespace Dune

#endif // DUNE_DPG_SYSTEM_ASSEMBLER_HH
