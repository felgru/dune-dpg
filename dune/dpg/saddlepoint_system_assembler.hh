// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_SADDLEPOINT_SYSTEM_ASSEMBLER_HH
#define DUNE_DPG_SADDLEPOINT_SYSTEM_ASSEMBLER_HH

#include <functional>
#include <list>
#include <map>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

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
#include <boost/fusion/sequence/intrinsic/value_at.hpp>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include "assemble_helper.hh"
#include "assemble_types.hh"
#include "bilinearform.hh"
#include "innerproduct.hh"
#include "quadrature.hh"
#include "linearform.hh"
#include "localevaluation.hh"

namespace Dune {

/**
 * \brief This constructs the matrix and vector of a DPG system.
 *
 * \tparam TSpaces        tuple of test spaces
 * \tparam SSpaces        tuple of solution spaces
 * \tparam BilinForm      bilinear form describing the system
 * \tparam InProduct      inner product of the test space
 */
template<class TSpaces, class SolSpaces,
         class BilinForm, class InProduct>
class SaddlepointSystemAssembler
{
public:
  typedef TSpaces TestSpaces;
  typedef SolSpaces SolutionSpaces;
  //! tuple type for the local views of the test spaces
  typedef typename boost::fusion::result_of::as_vector<
      typename boost::fusion::result_of::
      transform<TestSpaces, detail::getLocalView>::type>::type TestLocalViews;
  //! tuple type for the local views of the solution spaces
  typedef typename boost::fusion::result_of::as_vector<
      typename boost::fusion::result_of::
      transform<SolutionSpaces, detail::getLocalView>::type
      >::type SolutionLocalViews;
  //! type of the bilinear form describing this DPG system
  typedef BilinForm BilinearForm;
  //! type of the inner product on the test spaces
  typedef InProduct InnerProduct;

public:
  SaddlepointSystemAssembler () = delete;

  /**
   * \brief constructor for SaddlepointSystemAssembler
   *
   * \note For your convenience, use make_SaddlepointSystemAssembler()
   *       instead.
   */
  constexpr SaddlepointSystemAssembler (TestSpaces     testSpaces,
                                        SolutionSpaces solutionSpaces,
                                        BilinForm      bilinearForm,
                                        InProduct      innerProduct)
             : testSpaces(testSpaces),
               solutionSpaces(solutionSpaces),
               bilinearForm(detail::make_BilinearForm<SaddlepointFormulation>
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
   * \param[out] matrix        the matrix of the DPG system
   * \param[out] rhs           the rhs vector of the DPG system
   * \param[in]  rhsLinearForm the linear form describing the rhs
   */
  template <class LinearForm>
  void assembleSystem(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                      BlockVector<FieldVector<double,1> >& rhs,
                      LinearForm& rhsLinearForm);

  /**
   * \brief The same as assembleSystem but it only assembles the matrix.
   */
  void assembleMatrix(BCRSMatrix<FieldMatrix<double,1,1> >& matrix);

  /**
   * \brief The same as assembleSystem but it only assembles the rhs.
   */
  template <class LinearForm>
  void assembleRhs(BlockVector<FieldVector<double,1> >& rhs,
                   LinearForm& rhsLinearForm);

  /**
   * \brief Apply Dirichlet boundary values to a solution space
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
  void applyDirichletBoundary(
                              BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                              BlockVector<FieldVector<double,1> >& rhs,
                              const std::vector<bool>& dirichletNodes,
                              const ValueType& value);

  /**
   * \brief Apply Dirichlet boundary values to a solution space
   *
   * \param[in,out] matrix      the matrix of the DPG system
   * \param[in] dirichletNodes  true marks the dofs in the Dirichlet boundary
   * \param[in] value           the Dirichlet boundary value
   * \tparam spaceIndex  the index of the test space on which we apply
   *                     the boundary data
   * \tparam ValueType   we take either constants or functions for \p value
   */
  template <size_t spaceIndex, class ValueType>
  void applyDirichletBoundaryToMatrix(
                              BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                              const std::vector<bool>& dirichletNodes,
                              const ValueType& value);

  /**
   * \brief Apply Dirichlet boundary values to a solution space
   *
   * \param[in,out] rhs         the rhs vector of the DPG system
   * \param[in] dirichletNodes  true marks the dofs in the Dirichlet boundary
   * \param[in] value           the Dirichlet boundary value
   * \tparam spaceIndex  the index of the solution space on which we apply
   *                     the boundary data
   * \tparam ValueType   we take either constants or functions for \p value
   */
  template <size_t spaceIndex, class ValueType>
  void applyDirichletBoundaryToRhs(
                              BlockVector<FieldVector<double,1> >& rhs,
                              const std::vector<bool>& dirichletNodes,
                              const ValueType& value);

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

  template<class TestZip, class SolutionZip, class CP, class CPM>
  static inline void copy_local_matrix_to_global
      (const BilinearForm& bilinearForm,
       const InnerProduct& innerProduct,
       const TestZip& testZip,
       const SolutionZip& solutionZip,
       CP& cp,
       CPM& cpm);

  TestSpaces     testSpaces;
  SolutionSpaces solutionSpaces;
  BilinearForm   bilinearForm;
  InnerProduct   innerProduct;
};

/**
 * \brief Creates a SaddlepointSystemAssembler for a saddlepoint formulation,
 *        deducing the target type from the types of arguments.
 *
 * \param testSpaces     a tuple of test spaces
 * \param solutionSpaces a tuple of solution spaces
 * \param bilinearForm   the bilinear form describing the DPG system
 * \param innerProduct   the inner product of the test spaces
 */
template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct>
auto make_SaddlepointSystemAssembler(TestSpaces     testSpaces,
                                     SolutionSpaces solutionSpaces,
                                     BilinearForm   bilinearForm,
                                     InnerProduct   innerProduct)
    -> SaddlepointSystemAssembler<TestSpaces, SolutionSpaces,
                       BilinearForm, InnerProduct>
{
  return SaddlepointSystemAssembler<TestSpaces, SolutionSpaces,
                         BilinearForm, InnerProduct>
                      (testSpaces,
                       solutionSpaces,
                       bilinearForm,
                       innerProduct);
}


template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct>
template <class LinearForm>
void SaddlepointSystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct>::
assembleSystem(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
               BlockVector<FieldVector<double,1> >& rhs,
               LinearForm& rhsLinearForm)
{
  assembleMatrix(matrix);
  assembleRhs   (rhs, rhsLinearForm);
}

template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct>
void SaddlepointSystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct>::
assembleMatrix(BCRSMatrix<FieldMatrix<double,1,1> >& matrix)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  constexpr bool isSaddlepoint = true;

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  size_t globalTotalTestSize = 0;

  globalTotalTestSize =
      fold(zip(globalTestSpaceOffsets, testSpaces),
           (size_t)0, globalOffsetHelper());

  size_t globalTotalSolutionSize =
      fold(zip(globalSolutionSpaceOffsets, solutionSpaces),
           globalTotalTestSize, globalOffsetHelper());

  globalTotalSolutionSize -= globalSolutionSpaceOffsets[0];

  const auto n = globalTotalSolutionSize + globalTotalTestSize;

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  MatrixIndexSet occupationPattern;
  occupationPattern.resize(n, n);
  bilinearForm.template getOccupationPattern<isSaddlepoint>
               (occupationPattern,
                0, globalTotalTestSize);
  innerProduct.getOccupationPattern(occupationPattern);

  /* Add the diagonal of the matrix, since we need it for the
   * interpolation of boundary values. */
  for (size_t i=0; i<n; i++)
  {
    occupationPattern.add(i, i);
  }
  occupationPattern.exportIdx(matrix);

  // Set all entries to zero
  matrix = 0;

  // Views on the FE bases on a single element
  auto solutionLocalViews = as_vector(transform(solutionSpaces,
                                                getLocalView()));
  auto testLocalViews     = as_vector(transform(testSpaces,
                                                getLocalView()));

  auto solutionLocalIndexSets = as_vector(transform(solutionSpaces,
                                                    getLocalIndexSet()));
  auto testLocalIndexSets     = as_vector(transform(testSpaces,
                                                    getLocalIndexSet()));

  for(const auto& e : elements(gridView)) {

    for_each(solutionLocalViews, applyBind<decltype(e)>(e));
    for_each(testLocalViews, applyBind<decltype(e)>(e));

    for_each(zip(solutionLocalIndexSets, solutionLocalViews),
             make_fused_procedure(bindLocalIndexSet()));
    for_each(zip(testLocalIndexSets, testLocalViews),
             make_fused_procedure(bindLocalIndexSet()));

    bilinearForm.bind(testLocalViews, solutionLocalViews);
    innerProduct.bind(testLocalViews);

    // Now let's get the element stiffness matrix and the Gram matrix
    // for the test space.
    Matrix<FieldMatrix<double,1,1> > bfElementMatrix;
    Matrix<FieldMatrix<double,1,1> > ipElementMatrix;

    bilinearForm.getLocalMatrix(bfElementMatrix);
    innerProduct.getLocalMatrix(ipElementMatrix);


    // Add element stiffness matrix onto the global stiffness matrix
    auto cp = fused_procedure<localToGlobalCopier<decltype(ipElementMatrix),
                   typename std::remove_reference<decltype(matrix)>::type> >
                (localToGlobalCopier<decltype(ipElementMatrix),
                   typename std::remove_reference<decltype(matrix)>::type>
                                        (ipElementMatrix, matrix));
    auto cpm = fused_procedure<localToGlobalCopier<decltype(bfElementMatrix),
                    typename std::remove_reference<decltype(matrix)>::type,
                    isSaddlepoint> >
                 (localToGlobalCopier<decltype(bfElementMatrix),
                    typename std::remove_reference<decltype(matrix)>::type,
                    isSaddlepoint>
                                        (bfElementMatrix, matrix));

    /* copy every local submatrix indexed by a pair of indices from
     * bfIndices and ipIndices exactly once. */
    auto testZip = zip(testLocalViews,
                       testLocalIndexSets,
                       bilinearForm.getLocalTestSpaceOffsets(),
                       globalTestSpaceOffsets);
    auto solutionZip = zip(solutionLocalViews,
                           solutionLocalIndexSets,
                           bilinearForm.getLocalSolutionSpaceOffsets(),
                           globalSolutionSpaceOffsets);

    copy_local_matrix_to_global
      (bilinearForm,
       innerProduct,
       testZip,
       solutionZip,
       cp, cpm);
  }
}

template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct>
template <class LinearForm>
void SaddlepointSystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct>::
assembleRhs(BlockVector<FieldVector<double,1> >& rhs,
            LinearForm& rhsLinearForm)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  /* set up global offsets */
  size_t globalSpaceOffsets[std::tuple_size<TestSpaces>::value
                            + std::tuple_size<SolutionSpaces>::value];
  const size_t globalTotalSpaceSize =
      fold(zip(globalSpaceOffsets, join(testSpaces, solutionSpaces)),
           (size_t)0, globalOffsetHelper());

  // set rhs to correct length -- the total number of basis vectors in the bases
  rhs.resize(globalTotalSpaceSize);

  // Set all entries to zero
  rhs = 0;

  // Views on the FE bases on a single element
  auto localViews
        = as_vector(transform(join(testSpaces, solutionSpaces),
                              getLocalView()));

  auto localIndexSets
        = as_vector(transform(join(testSpaces, solutionSpaces),
                              getLocalIndexSet()));


  for(const auto& e : elements(gridView)) {

    for_each(localViews, applyBind<decltype(e)>(e));

    for_each(zip(localIndexSets, localViews),
             make_fused_procedure(bindLocalIndexSet()));

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> > localRhs;

    rhsLinearForm.bind(localViews);
    rhsLinearForm.getLocalVector(localRhs);

    auto cp = fused_procedure<localToGlobalRHSCopier<decltype(localRhs),
                   typename std::remove_reference<decltype(rhs)>::type> >
                (localToGlobalRHSCopier<decltype(localRhs),
                   typename std::remove_reference<decltype(rhs)>::type>
                                        (localRhs, rhs));

    /* copy every local subvector indexed by an index from
     * lfIndices exactly once. */
    auto testZip = zip(localViews,
                       localIndexSets,
                       rhsLinearForm.getLocalSpaceOffsets(),
                       globalSpaceOffsets);
    /* create set of indices to loop over when copying the local matrices
     * into the global one.
     */
    typedef typename boost::mpl::fold<
        typename boost::mpl::transform<
            /* This as_vector is probably not needed for boost::fusion 1.58
             * or higher. */
            typename result_of::as_vector<typename std::remove_reference<
                  decltype(rhsLinearForm.getTerms())>::type
                >::type
          , mpl::first<boost::mpl::_1>
          >::type
      , boost::mpl::set0<>
      , boost::mpl::insert<boost::mpl::_1,boost::mpl::_2>
      >::type LFIndices;

    auto lfIndices = LFIndices{};
    for_each(lfIndices,
             localToGlobalRHSCopyHelper<decltype(testZip),
                                        decltype(cp)>
                                       (testZip, cp));
  }
}


template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct>
template <size_t spaceIndex, class ValueType>
void SaddlepointSystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct>::
applyDirichletBoundary
                      (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                       BlockVector<FieldVector<double,1> >& rhs,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& boundaryValue)
{
  applyDirichletBoundaryToMatrix<spaceIndex>
          (matrix, dirichletNodes, boundaryValue);
  applyDirichletBoundaryToRhs<spaceIndex>
          (rhs,    dirichletNodes, boundaryValue);
}

template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct>
template <size_t spaceIndex, class ValueType>
void SaddlepointSystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct>::
applyDirichletBoundaryToMatrix
                      (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& boundaryValue)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  const size_t spaceSize =
        std::get<spaceIndex>(solutionSpaces).size();

  size_t globalOffset;
  {
    /* set up global offsets */
    size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
    size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
    size_t globalTotalTestSize = 0;

    globalTotalTestSize =
        fold(zip(globalTestSpaceOffsets, testSpaces),
             (size_t)0, globalOffsetHelper());

    fold(zip(globalSolutionSpaceOffsets, solutionSpaces),
         globalTotalTestSize, globalOffsetHelper());

    globalOffset = globalSolutionSpaceOffsets[spaceIndex];
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
         class BilinearForm, class InnerProduct>
template <size_t spaceIndex, class ValueType>
void SaddlepointSystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct>::
applyDirichletBoundaryToRhs
                      (BlockVector<FieldVector<double,1> >& rhs,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& boundaryValue)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  const size_t spaceSize =
        std::get<spaceIndex>(solutionSpaces).size();

  size_t globalOffset;
  {
    /* set up global offsets */
    size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
    size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
    size_t globalTotalTestSize = 0;

    globalTotalTestSize =
        fold(zip(globalTestSpaceOffsets, testSpaces),
             (size_t)0, globalOffsetHelper());

    fold(zip(globalSolutionSpaceOffsets, solutionSpaces),
         globalTotalTestSize, globalOffsetHelper());

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

}


template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct>
template<size_t spaceIndex, unsigned int dim>
void SaddlepointSystemAssembler<TestSpaces, SolutionSpaces,
                     BilinearForm, InnerProduct>::
defineCharacteristicFaces(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                          BlockVector<FieldVector<double,1> >& rhs,
                          const FieldVector<double,dim>& beta,
                          double delta)
{
  static_assert(std::is_same<BilinearForm, TestSpaces>::value,
                "defineCharacteristicFaces not implemented "
                "for Saddlepointformulation.");
}


template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct>
template <size_t spaceIndex, class MinInnerProduct, unsigned int dim>
void SaddlepointSystemAssembler<TestSpaces, SolutionSpaces,
                                BilinearForm, InnerProduct>::
applyMinimization
            (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
             MinInnerProduct minInnerProduct,
             FieldVector<double, dim> beta,
             double delta,
             double epsilon)
{
  static_assert(std::is_same<BilinearForm, TestSpaces>::value,
                "applyMinimization not implemented "
                "for Saddlepointformulation ");
}

template<class TestSpaces, class SolutionSpaces,
         class BilinearForm, class InnerProduct>
template<class TestZip, class SolutionZip, class CP, class CPM>
inline void SaddlepointSystemAssembler<TestSpaces, SolutionSpaces,
                                       BilinearForm, InnerProduct>::
copy_local_matrix_to_global
    (const BilinearForm& bilinearForm,
     const InnerProduct& innerProduct,
     const TestZip& testZip,
     const SolutionZip& solutionZip,
     CP& cp,
     CPM& cpm)
{
  using namespace boost::fusion;
  using namespace detail;

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

  for_each(ipIndices,
           localToGlobalCopyHelper<TestZip,
                                   TestZip,
                                   CP>
                                  (testZip, testZip, cp));
  for_each(bfIndices,
           localToGlobalCopyHelper<SolutionZip,
                                   TestZip,
                                   CPM>
                                  (solutionZip, testZip, cpm));
}


} // end namespace Dune

#endif // DUNE_DPG_SADDLEPOINT_SYSTEM_ASSEMBLER_HH
