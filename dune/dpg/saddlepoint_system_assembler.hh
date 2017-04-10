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

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>
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
 * \tparam BilinForm      bilinear form describing the system
 * \tparam InProduct      inner product of the test space
 */
template<class BilinForm, class InProduct>
class SaddlepointSystemAssembler
{
public:
  typedef typename BilinForm::TestSpaces TestSpaces;
  typedef typename BilinForm::SolutionSpaces SolutionSpaces;
  //! tuple type for the local views of the test spaces
  typedef detail::getLocalViews_t<TestSpaces>  TestLocalViews;
  //! tuple type for the local views of the solution spaces
  typedef detail::getLocalViews_t<SolutionSpaces>  SolutionLocalViews;
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
  constexpr SaddlepointSystemAssembler (const BilinForm&   bilinearForm,
                                        const InProduct&   innerProduct)
             : testSpaces(bilinearForm.getTestSpaces()),
               solutionSpaces(bilinearForm.getSolutionSpaces()),
               bilinearForm(bilinearForm),
               innerProduct(innerProduct)
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
  std::shared_ptr<TestSpaces> getTestSpaces() const
  { return testSpaces; }

  /**
   * \brief Does exactly what it says on the tin.
   */
  std::shared_ptr<SolutionSpaces> getSolutionSpaces() const
  { return solutionSpaces; }

private:

  /* create sets of index pairs to loop over.
   * This will be used later, when copying the local matrices into
   * the global one.
   */
  static constexpr auto bfIndices() {
    namespace hana = boost::hana;
    using BilinearTerms = typename BilinearForm::Terms;
    auto bfIndices = hana::to<hana::set_tag>(
        hana::transform(hana::to<hana::tuple_tag>(
            hana::make_range(hana::int_c<0>,
                hana::int_c<std::tuple_size<BilinearTerms>::value>)),
          [](auto i) -> auto {
            using Term = std::tuple_element_t<i.value, BilinearTerms>;
            return hana::tuple<std::tuple_element_t<0, Term>,
                               std::tuple_element_t<1, Term>>{};
          }));
    return hana::to<hana::tuple_tag>(bfIndices);
  }

  using BFIndices = decltype(bfIndices());

  static constexpr auto ipIndices() {
    namespace hana = boost::hana;
    using InnerProductTerms = typename InnerProduct::Terms;
    auto ipIndices = hana::to<hana::set_tag>(
        hana::transform(hana::to<hana::tuple_tag>(
            hana::make_range(hana::int_c<0>,
                hana::int_c<std::tuple_size<InnerProductTerms>::value>)),
          [](auto i) -> auto {
            using Term = std::tuple_element_t<i.value, InnerProductTerms>;
            return hana::tuple<std::tuple_element_t<0, Term>,
                               std::tuple_element_t<1, Term>>{};
          }));
    return hana::to<hana::tuple_tag>(ipIndices);
  }

  // TODO: IPIndices seems to have some type_c too much.
  using IPIndices = decltype(ipIndices());

  std::shared_ptr<TestSpaces>     testSpaces;
  std::shared_ptr<SolutionSpaces> solutionSpaces;
  BilinearForm                    bilinearForm;
  InnerProduct                    innerProduct;
};

/**
 * \brief Creates a SaddlepointSystemAssembler for a saddlepoint formulation,
 *        deducing the target type from the types of arguments.
 *
 * \param bilinearForm   the bilinear form describing the DPG system
 * \param innerProduct   the inner product of the test spaces
 */
template<class BilinearForm, class InnerProduct>
auto make_SaddlepointSystemAssembler(const BilinearForm&  bilinearForm,
                                     const InnerProduct&  innerProduct)
    -> SaddlepointSystemAssembler<BilinearForm, InnerProduct>
{
  return SaddlepointSystemAssembler<BilinearForm, InnerProduct>
                      (bilinearForm,
                       innerProduct);
}


template<class BilinearForm, class InnerProduct>
template <class LinearForm>
void SaddlepointSystemAssembler<BilinearForm, InnerProduct>::
assembleSystem(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
               BlockVector<FieldVector<double,1> >& rhs,
               LinearForm& rhsLinearForm)
{
  assembleMatrix(matrix);
  assembleRhs   (rhs, rhsLinearForm);
}

template<class BilinearForm, class InnerProduct>
void SaddlepointSystemAssembler<BilinearForm, InnerProduct>::
assembleMatrix(BCRSMatrix<FieldMatrix<double,1,1> >& matrix)
{
  using namespace Dune::detail;

  constexpr bool isSaddlepoint = true;

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(*testSpaces).gridView();

  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  const size_t globalTotalTestSize = computeOffsets(globalTestSpaceOffsets,
                                                    *testSpaces);

  const size_t globalTotalSolutionSize =
      computeOffsets(globalSolutionSpaceOffsets, *solutionSpaces,
                     globalTotalTestSize) - globalTotalTestSize;

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
  auto solutionLocalViews = getLocalViews(*solutionSpaces);
  auto testLocalViews     = getLocalViews(*testSpaces);

  auto solutionLocalIndexSets = getLocalIndexSets(*solutionSpaces);
  auto testLocalIndexSets = getLocalIndexSets(*testSpaces);

  for(const auto& e : elements(gridView)) {

    bindLocalViews(solutionLocalViews, e);
    bindLocalViews(testLocalViews, e);

    bindLocalIndexSets(solutionLocalIndexSets, solutionLocalViews);
    bindLocalIndexSets(testLocalIndexSets, testLocalViews);

    bilinearForm.bind(testLocalViews, solutionLocalViews);
    innerProduct.bind(testLocalViews);

    // Now let's get the element stiffness matrix and the Gram matrix
    // for the test space.
    Matrix<FieldMatrix<double,1,1> > bfElementMatrix;
    Matrix<FieldMatrix<double,1,1> > ipElementMatrix;

    bilinearForm.getLocalMatrix(bfElementMatrix);
    innerProduct.getLocalMatrix(ipElementMatrix);

    const auto& localTestSpaceOffsets
        = bilinearForm.getLocalTestSpaceOffsets();
    const auto& localSolutionSpaceOffsets
        = bilinearForm.getLocalSolutionSpaceOffsets();

    // Add element stiffness matrix onto the global stiffness matrix
    /* copy every local submatrix indexed by a pair of indices from
     * bfIndices and ipIndices exactly once. */
    copyLocalToGlobalMatrix<IPIndices>(
        ipElementMatrix,
        matrix,
        testLocalIndexSets,
        localTestSpaceOffsets,
        globalTestSpaceOffsets,
        testLocalIndexSets,
        localTestSpaceOffsets,
        globalTestSpaceOffsets);
    copyLocalToGlobalMatrixSymmetric<BFIndices>(
        bfElementMatrix,
        matrix,
        testLocalIndexSets,
        localTestSpaceOffsets,
        globalTestSpaceOffsets,
        solutionLocalIndexSets,
        localSolutionSpaceOffsets,
        globalSolutionSpaceOffsets);
  }
}

template<class BilinearForm, class InnerProduct>
template <class LinearForm>
void SaddlepointSystemAssembler<BilinearForm, InnerProduct>::
assembleRhs(BlockVector<FieldVector<double,1> >& rhs,
            LinearForm& rhsLinearForm)
{
  using namespace Dune::detail;

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(*testSpaces).gridView();

  /* set up global offsets */
  size_t globalSpaceOffsets[std::tuple_size<TestSpaces>::value
                            + std::tuple_size<SolutionSpaces>::value];
  const size_t globalTotalSpaceSize =
      computeOffsets(globalSpaceOffsets,
                     std::tuple_cat(*testSpaces, *solutionSpaces));

  // set rhs to correct length -- the total number of basis vectors in the bases
  rhs.resize(globalTotalSpaceSize);

  // Set all entries to zero
  rhs = 0;

  // Views on the FE bases on a single element
  auto localViews
        = getLocalViews(std::tuple_cat(*testSpaces, *solutionSpaces));

  auto localIndexSets
        = getLocalIndexSets(std::tuple_cat(*testSpaces, *solutionSpaces));


  for(const auto& e : elements(gridView)) {

    bindLocalViews(localViews, e);
    bindLocalIndexSets(localIndexSets, localViews);

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> > localRhs;

    rhsLinearForm.bind(localViews);
    rhsLinearForm.getLocalVector(localRhs);

    /* create set of indices to loop over when copying the local matrices
     * into the global one.
     */
    namespace hana = boost::hana;
    using LinearTerms = std::decay_t<decltype(rhsLinearForm.getTerms())>;
    auto lfIndices = hana::to<hana::set_tag>(
        hana::transform(hana::to<hana::tuple_tag>(
            hana::make_range(hana::int_c<0>,
                hana::int_c<std::tuple_size<LinearTerms>::value>)),
          [](auto i) -> auto {
            using Term = std::tuple_element_t<i.value, LinearTerms>;
            return hana::type_c<std::tuple_element_t<0, Term>>;
          }));
    using LFIndices = decltype(hana::to<hana::tuple_tag>(lfIndices));

    // Add local right-hand side onto the global right-hand side
    /* copy every local subvector indexed by an index from
     * lfIndices exactly once. */
    copyLocalToGlobalVector<LFIndices>(
        localRhs,
        rhs,
        localIndexSets,
        rhsLinearForm.getLocalSpaceOffsets(),
        globalSpaceOffsets);
  }
}


template<class BilinearForm, class InnerProduct>
template <size_t spaceIndex, class ValueType>
void SaddlepointSystemAssembler<BilinearForm, InnerProduct>::
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

template<class BilinearForm, class InnerProduct>
template <size_t spaceIndex, class ValueType>
void SaddlepointSystemAssembler<BilinearForm, InnerProduct>::
applyDirichletBoundaryToMatrix
                      (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& boundaryValue)
{
  const size_t spaceSize =
        std::get<spaceIndex>(*solutionSpaces).size();

  const size_t globalOffset
      = detail::computeOffset<spaceIndex>(*solutionSpaces,
              detail::computeOffset<std::tuple_size<TestSpaces>::value>
                                   (*testSpaces));

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

template<class BilinearForm, class InnerProduct>
template <size_t spaceIndex, class ValueType>
void SaddlepointSystemAssembler<BilinearForm, InnerProduct>::
applyDirichletBoundaryToRhs
                      (BlockVector<FieldVector<double,1> >& rhs,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& boundaryValue)
{
  const size_t spaceSize =
        std::get<spaceIndex>(*solutionSpaces).size();

  const size_t globalOffset
      = detail::computeOffset<spaceIndex>(*solutionSpaces,
              detail::computeOffset<std::tuple_size<TestSpaces>::value>
                                   (*testSpaces));

  // Set Dirichlet values
  for (size_t i=0; i<spaceSize; i++)
  {
    if (dirichletNodes[i])
    {
      /* TODO: Needs adaptation when value is a function. */
      rhs[globalOffset+i] = detail::evaluateFactor(boundaryValue, i);
    }
  }

}


template<class BilinearForm, class InnerProduct>
template<size_t spaceIndex, unsigned int dim>
void SaddlepointSystemAssembler<BilinearForm, InnerProduct>::
defineCharacteristicFaces(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                          BlockVector<FieldVector<double,1> >& rhs,
                          const FieldVector<double,dim>& beta,
                          double delta)
{
  static_assert(std::is_same<BilinearForm, TestSpaces>::value,
                "defineCharacteristicFaces not implemented "
                "for Saddlepointformulation.");
}


template<class BilinearForm, class InnerProduct>
template <size_t spaceIndex, class MinInnerProduct, unsigned int dim>
void SaddlepointSystemAssembler<BilinearForm, InnerProduct>::
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


} // end namespace Dune

#endif // DUNE_DPG_SADDLEPOINT_SYSTEM_ASSEMBLER_HH
