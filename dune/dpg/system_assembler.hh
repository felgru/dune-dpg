// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_SYSTEM_ASSEMBLER_HH
#define DUNE_DPG_SYSTEM_ASSEMBLER_HH

#include <functional>
#include <list>
#include <map>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

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
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/container/set/convert.hpp>
#include <boost/fusion/algorithm/auxiliary/copy.hpp>
#include <boost/fusion/algorithm/transformation/join.hpp>
#include <boost/fusion/algorithm/transformation/transform.hpp>
#include <boost/fusion/algorithm/transformation/zip.hpp>
#include <boost/fusion/algorithm/iteration/accumulate.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/functional/generation/make_fused_procedure.hpp>

#include <dune/geometry/quadraturerules.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
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
template <class LocalViewTest, class LocalVolumeTerm, class TestSpace>
void getVolumeTerm(const LocalViewTest& localViewTest,
                   BlockVector<FieldVector<double,1> >& localRhs,
                   LocalVolumeTerm&& localVolumeTerm,
                   TestSpace&)
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

  /* TODO: Quadrature order is only good enough for a constant localVolumeTerm. */
  int quadratureOrder = localFiniteElementTest.localBasis().order();

  constexpr bool useSubsampledQuadrature =
      is_SubsampledFiniteElement<TestSpace>::value;

  if(useSubsampledQuadrature) {
    constexpr int s = numberOfSamples<TestSpace>::value;

    const QuadratureRule<double, dim>& quadSection =
          QuadratureRules<double, dim>::rule(element.type(), quadratureOrder);
    const SubsampledQuadratureRule<double, s, dim> quad(quadSection);

    for ( size_t pt=0; pt < quad.size(); pt++ ) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();

      // The multiplicative factor in the integral transformation formula
      const double integrationElement =
              element.geometry().integrationElement(quadPos);

      double functionValue = localVolumeTerm(quadPos);

      std::vector<FieldVector<double,1> > shapeFunctionValues;
      localFiniteElementTest.localBasis().evaluateFunction(quadPos,
                                                           shapeFunctionValues);

      for (size_t i=0; i<localRhs.size(); i++)
        localRhs[i] += shapeFunctionValues[i] * functionValue
                     * quad[pt].weight() * integrationElement;

    }
  } else {
    const QuadratureRule<double, dim>& quad =
        QuadratureRules<double, dim>::rule(element.type(), quadratureOrder);

    for ( size_t pt=0; pt < quad.size(); pt++ ) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();

      // The multiplicative factor in the integral transformation formula
      const double integrationElement =
              element.geometry().integrationElement(quadPos);

      double functionValue = localVolumeTerm(quadPos);

      std::vector<FieldVector<double,1> > shapeFunctionValues;
      localFiniteElementTest.localBasis().evaluateFunction(quadPos,
                                                           shapeFunctionValues);

      for (size_t i=0; i<localRhs.size(); i++)
        localRhs[i] += shapeFunctionValues[i] * functionValue
                     * quad[pt].weight() * integrationElement;

    }
  }

}

namespace detail {
struct getVolumeTermHelper
{
  template<class Seq>
  void operator()(const Seq& seq) const
  {
    using namespace boost::fusion;

    // TODO: this can probably be done more elegantly by sequence fusion.
    getVolumeTerm(at_c<0>(seq), at_c<1>(seq), at_c<2>(seq), at_c<3>(seq));
  }
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
               bilinearForm(detail::make_BilinearForm<FormulationType>
                                             (testSpaces,
                                              solutionSpaces,
                                              bilinearForm.getTerms())),
               innerProduct(make_InnerProduct(testSpaces,
                                              innerProduct.getTerms()))
  { }

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
                              double mu);

  /**
   * \brief Linearly interpolate dofs on characteristic faces.
   *
   * The inner degrees of freedom on the characteristic faces are
   * undefined in the DPG formulation.  This function assigns them
   * a well-defined value by interpolating between the end points of
   * the face.
   * This interpolation makes sure that we get a well-posed linear
   * system.
   *
   * \param[in,out] matrix  the matrix of the DPG system
   * \param[in,out] rhs     the rhs vector of the DPG system
   * \param[in]     beta    the transport direction
   * \param         delta   tolerance for numeric definition of
   *                           characteristic face
   */
  template <size_t spaceIndex, unsigned int dim>
  void defineCharacteristicFaces(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                                 BlockVector<FieldVector<double,1> >& rhs,
                                 const FieldVector<double,dim>& beta,
                                 double delta = 10e-10);

  template <size_t spaceIndex, class MinInnerProduct, unsigned int dim>
  void applyMinimization(
                      BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                      MinInnerProduct minInnerProduct,
                      FieldVector<double, dim> beta,
                      double delta = 10e-10,
                      double epsilon = 0);

  /**
   * \brief Does exactly what it says on the tin.
   */
  const TestSpaces& getTestSpaces() const
  { return testSpaces; }

  /**
   * \brief Does exactly what it says on the tin.
   */
  const SolutionSpaces& getSolutionSpaces() const
  { return solutionSpaces; }

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

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  /* TODO: make this a vector of pointers to localVolumeTerm */
  auto localVolumeTerms =
      as_vector(transform(volumeTerms,
                          getLocalVolumeTerm<GridView>(gridView)));


  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  size_t globalTotalTestSize = 0;

  if(isSaddlepoint) {
    globalTotalTestSize =
        fold(zip(globalTestSpaceOffsets, testSpaces),
             (size_t)0, globalOffsetHelper());
  } else { /* DPG formulation */
    for(size_t i=0; i<std::tuple_size<TestSpaces>::value; ++i)
    {
      globalTestSpaceOffsets[i] = 0;
    }
  }

  size_t globalTotalSolutionSize =
      fold(zip(globalSolutionSpaceOffsets, solutionSpaces),
           isSaddlepoint?globalTotalTestSize:0, globalOffsetHelper());

  globalTotalSolutionSize -= globalSolutionSpaceOffsets[0];

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
  occupationPattern.exportIdx(matrix);

  // set rhs to correct length -- the total number of basis vectors in the basis
  rhs.resize((isSaddlepoint?globalTotalTestSize:0) + globalTotalSolutionSize);

  // Set all entries to zero
  matrix = 0;
  rhs = 0;

  // Views on the FE bases on a single element
  auto solutionLocalView = as_vector(transform(solutionSpaces, getLocalView()));
  auto testLocalView     = as_vector(transform(testSpaces, getLocalView()));

  auto solutionLocalIndexSet = as_vector(transform(solutionSpaces,
                                                   getLocalIndexSet()));
  auto testLocalIndexSet     = as_vector(transform(testSpaces,
                                                   getLocalIndexSet()));

  /* create sets of index pairs to loop over.
   * This will be used later, when copying the local matrices into
   * the global one.
   */
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

  for(const auto& e : elements(gridView)) {

    for_each(solutionLocalView, applyBind<decltype(e)>(e));
    for_each(testLocalView, applyBind<decltype(e)>(e));

    for_each(zip(solutionLocalIndexSet, solutionLocalView),
             make_fused_procedure(bindLocalIndexSet()));
    for_each(zip(testLocalIndexSet, testLocalView),
             make_fused_procedure(bindLocalIndexSet()));

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

    /* copy every local submatrix indexed by a pair of indices from
     * bfIndices and ipIndices exactly once. */
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
                                decltype(localVolumeTerms)&,
                                decltype(testSpaces)&>;
    for_each(zip_view<RHSZipHelper>(RHSZipHelper(testLocalView,
                                                 localRhs,
                                                 localVolumeTerms,
                                                 testSpaces)),
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
                       const ValueType& boundaryValue,
                       const Spaces& spaces)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  constexpr bool isSaddlepoint =
        std::is_same<
             typename std::decay<FormulationType>::type
           , SaddlepointFormulation
        >::value;

  const size_t spaceSize =
        std::get<spaceIndex>(spaces).size();

  size_t globalOffset;
  {
    /* set up global offsets */
    size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
    size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
    size_t globalTotalTestSize = 0;

    if(isSaddlepoint) {
      globalTotalTestSize =
          fold(zip(globalTestSpaceOffsets, testSpaces),
               (size_t)0, globalOffsetHelper());
    } else { /* DPG formulation */
      for(size_t i=0; i<std::tuple_size<TestSpaces>::value; ++i)
      {
        globalTestSpaceOffsets[i] = 0;
      }
    }

    fold(zip(globalSolutionSpaceOffsets, solutionSpaces),
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
      rhs[globalOffset+i] = detail::evaluateFactor(boundaryValue,i);
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
          /* Zero out row and column to keep symmetry. */
          *cIt = 0;
          matrix[cIt.index()][globalOffset+i]=0;
        }
      }
    }

  }
}

template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct,
         class FormulationType>
template <size_t spaceIndex, unsigned int dim>
void SystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct, FormulationType>::
applyWeakBoundaryCondition
                    (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                     FieldVector<double, dim> beta,
                     double mu)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  static_assert(std::is_same<
             typename std::decay<FormulationType>::type
           , DPGFormulation
                            >::value,
                "applyWeakBoundaryConditions not implemented "
                "for Saddlepointformulation ");

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  fold(zip(globalSolutionSpaceOffsets, solutionSpaces),
       (size_t)0, globalOffsetHelper());
  size_t globalOffset = globalSolutionSpaceOffsets[spaceIndex];

  auto localView = std::get<spaceIndex>(solutionSpaces).localView();
  typedef decltype(localView) LocalView;
  auto localIndexSet = at_c<spaceIndex>(solutionSpaces).localIndexSet();

  for(const auto& e : elements(gridView))
  {
    localView.bind(e);
    localIndexSet.bind(localView);

    const auto& localFiniteElement = localView.tree().finiteElement();

    int quadratureOrder = 2*localFiniteElement.localBasis().order();

    size_t n = localFiniteElement.localBasis().size();

    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    // Set all matrix entries to zero
    elementMatrix.setSize(n,n);
    elementMatrix = 0;

    for (auto&& intersection : intersections(gridView, e))
    {
      if (intersection.boundary()) //if the intersection is at the (physical) boundary of the domain
      {
        const FieldVector<double,dim>& centerOuterNormal =
                intersection.centerUnitOuterNormal();

        if ((beta*centerOuterNormal) > -1e-10) //everywhere except inflow boundary
        {
          const QuadratureRule<double, dim-1>& quadFace =
                  QuadratureRules<double, dim-1>::rule(intersection.type(),
                                                       quadratureOrder);

          for (size_t pt=0; pt < quadFace.size(); pt++)
          {
            // position of the current quadrature point in the reference element (face!)
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
                        += (mu * solutionValues[i] * solutionValues[j])
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
}


template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct,
         class FormulationType>
template<size_t spaceIndex, unsigned int dim>
void SystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct, FormulationType>::
defineCharacteristicFaces(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                          BlockVector<FieldVector<double,1> >& rhs,
                          const FieldVector<double,dim>& beta,
                          double delta)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  static_assert(std::is_same<
             typename std::decay<FormulationType>::type
           , DPGFormulation
                            >::value,
                "defineCharacteristicFaces not implemented "
                "for Saddlepointformulation.");

  typedef typename std::tuple_element<spaceIndex,SolutionSpaces>::type::GridView GridView;
  GridView gridView = std::get<spaceIndex>(solutionSpaces).gridView();

  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  fold(zip(globalSolutionSpaceOffsets, solutionSpaces),
       (size_t)0, globalOffsetHelper());
  size_t globalOffset = globalSolutionSpaceOffsets[spaceIndex];

  size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

  auto solutionLocalView = at_c<spaceIndex>(solutionSpaces).localView();
  auto localIndexSet = at_c<spaceIndex>(solutionSpaces).localIndexSet();

  for(const auto& e : elements(gridView))
  {
    solutionLocalView.bind(e);
    localIndexSet.bind(solutionLocalView);

    const auto& localFiniteElement = solutionLocalView.tree().finiteElement();

    size_t n = localFiniteElement.localBasis().size();

    std::vector<bool>  characteristicFaces(e.subEntities(1), false);
    bool characteristicFound = false;

    for (auto&& intersection : intersections(gridView, e))
    {
      const bool characteristic =
          fabs(beta * intersection.centerUnitOuterNormal()) < delta;
      characteristicFaces[intersection.indexInInside()] = characteristic;
      characteristicFound = characteristicFound || characteristic;
    }

    if(characteristicFound)
    {
      std::map<size_t,std::list<std::pair<size_t,size_t>>> characteristicDOFs;
      std::vector<size_t> vertexDOFs(e.subEntities(dim));

      for (unsigned int i=0; i<n; i++)
      {
        if (localFiniteElement.localCoefficients().localKey(i).codim()==1)
            // edge DOFs
        {
          const size_t face =
              localFiniteElement.localCoefficients().localKey(i).subEntity();
          const size_t localIndex =
              localFiniteElement.localCoefficients().localKey(i).index();
          if(characteristicFaces[face])
            characteristicDOFs[face].emplace_back(i,localIndex);
        } else if (localFiniteElement.localCoefficients().localKey(i).codim()
                   == dim)
        {
          const size_t vertex =
              localFiniteElement.localCoefficients().localKey(i).subEntity();
          vertexDOFs[vertex] = i;
        }
        // Vertex DOFs are never characteristic because the corresponding
        // basis functions have support on at least two edges which can
        // never be both (almost) characteristic.
      }

      std::vector<std::pair<size_t,size_t>> endpoints(e.subEntities(1));
      /* TODO: this code only works on quadrilateral elements: */
      endpoints[0] = std::make_pair(vertexDOFs[0], vertexDOFs[2]);
      endpoints[1] = std::make_pair(vertexDOFs[1], vertexDOFs[3]);
      endpoints[2] = std::make_pair(vertexDOFs[0], vertexDOFs[1]);
      endpoints[3] = std::make_pair(vertexDOFs[2], vertexDOFs[3]);

      for (auto&& faceAndDOFs: characteristicDOFs)
      {
        size_t face;
        std::list<std::pair<size_t,size_t>> dofs;
        std::tie(face, dofs) = faceAndDOFs;
          size_t left, right;
          std::tie(left, right) = endpoints[face];
        for(auto&& dof: dofs)
        {
          auto row = localIndexSet.index(dof.first)[0];
          auto col = row;
          const size_t k = dofs.size()+1;

          /* replace the row of dof on characteristic face
           * by an interpolation of the two endpoints of the
           * characteristic face. */
          matrix[row+globalOffset][col+globalOffset] = -1;
          col = localIndexSet.index(left)[0];
          matrix[row+globalOffset][col+globalOffset]
              = (double)(k-dof.second-1)/k;
          col = localIndexSet.index(right)[0];
          matrix[row+globalOffset][col+globalOffset]
              = (double)(dof.second+1)/k;

          rhs[row+globalOffset] = 0;
        }
      }
    }
  }
}


template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct,
         class FormulationType>
template <size_t spaceIndex, class MinInnerProduct, unsigned int dim>
void SystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct, FormulationType>::
applyMinimization
            (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
             MinInnerProduct minInnerProduct,
             FieldVector<double, dim> beta,
             double delta,
             double epsilon)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  static_assert(std::is_same<
             typename std::decay<FormulationType>::type
           , DPGFormulation
                            >::value,
                "applyMinimization not implemented "
                "for Saddlepointformulation ");

  typedef typename std::tuple_element<spaceIndex,SolutionSpaces>::type::GridView GridView;
  GridView gridView = std::get<spaceIndex>(solutionSpaces).gridView();

  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  fold(zip(globalSolutionSpaceOffsets, solutionSpaces),
       (size_t)0, globalOffsetHelper());
  size_t globalOffset = globalSolutionSpaceOffsets[spaceIndex];

  size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

  // get local view for solution space (necessary if we want to use inner product) /TODO inefficient (why?)
  auto solutionLocalView = as_vector(transform(solutionSpaces, getLocalView()));

  auto localIndexSet = at_c<spaceIndex>(solutionSpaces).localIndexSet();

  bool epsilonSmallerDelta(epsilon<delta);

  for(const auto& e : elements(gridView))
  {
    for_each(solutionLocalView, applyBind<decltype(e)>(e));
    localIndexSet.bind(at_c<spaceIndex>(solutionLocalView));

    /* set up local offsets */
    fold(zip(localSolutionSpaceOffsets, solutionLocalView),
         (size_t)0, offsetHelper());

    const auto& localFiniteElement = at_c<spaceIndex>(solutionLocalView).tree().finiteElement();

    size_t n = localFiniteElement.localBasis().size();

    Matrix<FieldMatrix<double,1,1> > elementMatrix;

    minInnerProduct.bind(solutionLocalView);
    minInnerProduct.getLocalMatrix(elementMatrix);

    std::vector<bool> relevantFaces(e.subEntities(1), false);

    for (auto&& intersection : intersections(gridView, e))
    {
      const FieldVector<double,dim>& centerOuterNormal =
              intersection.centerUnitOuterNormal();
      //set relevant faces for almost characteristic faces
      relevantFaces[intersection.indexInInside()] = (std::abs(beta*centerOuterNormal) < delta);
    }

    std::vector<bool> relevantDOFs(n, false);

    for (unsigned int i=0; i<n; i++)
    {
      if (localFiniteElement.localCoefficients().localKey(i).codim()==0) // interior DOFs
      {
        relevantDOFs[i] = true;
      }
      else if (localFiniteElement.localCoefficients().localKey(i).codim()==1 and epsilonSmallerDelta) // edge DOFs
      {
        relevantDOFs[i] = relevantFaces[localFiniteElement.localCoefficients().localKey(i).subEntity()];
      }
      // Vertex DOFs are never relevant because the corresponding
      // basis functions have support on at least two edges which can
      // never be both (almost) characteristic.
    }

    for (size_t i=0; i<n; i++)
    {
      if (relevantDOFs[i])
      {
        auto row = localIndexSet.index(i)[0];
        for (size_t j=0; j<n; j++)
        {
          auto col = localIndexSet.index(j)[0];
          matrix[row+globalOffset][col+globalOffset]
                       += elementMatrix[localSolutionSpaceOffsets[spaceIndex]+i]
                                       [localSolutionSpaceOffsets[spaceIndex]+j];
        }
      }
    }
  }
}







} // end namespace Dune

#endif // DUNE_DPG_SYSTEM_ASSEMBLER_HH
