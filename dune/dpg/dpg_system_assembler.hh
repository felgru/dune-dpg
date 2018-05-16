// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_SYSTEM_ASSEMBLER_DPG_HH
#define DUNE_DPG_SYSTEM_ASSEMBLER_DPG_HH

#include <iterator>
#include <list>
#include <map>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

#include <boost/hana.hpp>

#include <dune/common/version.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dune/functions/functionspacebases/interpolate.hh>

#include <dune/dpg/functions/localindexsetiteration.hh>
#include "assemble_helper.hh"
#include "assemble_types.hh"
#include "bilinearform.hh"
#include "innerproduct.hh"
#include "quadrature.hh"
#include "linearform.hh"
#include "localevaluation.hh"
#include "testspace_coefficient_matrix.hh"


namespace Dune {

namespace detail {

struct Unbuffered {
  template<class BilinearForm, class InnerProduct>
  using TestspaceCoefficientMatrix
    = Dune::UnbufferedTestspaceCoefficientMatrix<BilinearForm, InnerProduct>;
};

struct Buffered {
  template<class BilinearForm, class InnerProduct>
  using TestspaceCoefficientMatrix
    = Dune::BufferedTestspaceCoefficientMatrix<BilinearForm, InnerProduct>;
};
}

/**
 * \brief This constructs the system matrix and righthandside vector of a DPG system.
 *
 * \tparam BilinForm         bilinear form
 * \tparam InnProduct        inner product
 * \tparam BufferPolicy      Buffered or Unbuffered
 */
template<class BilinForm, class InnProduct, class BufferPolicy>
class DPGSystemAssembler
{
  static_assert(std::is_same<BufferPolicy, detail::Unbuffered>::value
      or (uses_only_constant_coefficients_v<typename BilinForm::Terms>
          and uses_only_constant_coefficients_v<typename InnProduct::Terms>),
      "The GeometryBuffer in DPGSystemAssembler only works when all"
      " coefficients in the bilinear form and inner product are"
      " constant.");
public:
  using BilinearForm = BilinForm;
  using InnerProduct = InnProduct;
  using TestspaceCoefficientMatrix
    = typename BufferPolicy::template
        TestspaceCoefficientMatrix<BilinearForm, InnerProduct>;
  using TestSearchSpacesPtr = typename BilinearForm::TestSpacesPtr;
  using SolutionSpacesPtr = typename BilinearForm::SolutionSpacesPtr;
  using TestSearchSpaces = typename BilinearForm::TestSpaces;
  using SolutionSpaces = typename BilinearForm::SolutionSpaces;


  //! tuple type for the local views of the test spaces
  using TestLocalViews = detail::getLocalViews_t<TestSearchSpaces>;
  //! tuple type for the local views of the solution spaces
  using SolutionLocalViews = detail::getLocalViews_t<SolutionSpaces>;

  DPGSystemAssembler () = delete;
  /**
   * \brief constructor for DPGSystemAssembler
   *
   * \note For your convenience, use make_DPGSystemAssembler()
   *       instead.
   */
  constexpr DPGSystemAssembler (BilinearForm&  bilinearForm,
                                InnerProduct&  innerProduct)
             : testSearchSpaces_(bilinearForm.getTestSpaces()),
               solutionSpaces_(bilinearForm.getSolutionSpaces()),
               bilinearForm_(bilinearForm),
               testspaceCoefficientMatrix_(bilinearForm,
                                           innerProduct)
  { }

  template <class GeometryBuffer>
  constexpr DPGSystemAssembler (BilinearForm&      bilinearForm,
                                InnerProduct&      innerProduct,
                                GeometryBuffer&    geometryBuffer
                               )
             : testSearchSpaces_(bilinearForm.getTestSpaces()),
               solutionSpaces_(bilinearForm.getSolutionSpaces()),
               bilinearForm_(bilinearForm),
               testspaceCoefficientMatrix_(bilinearForm,
                                           innerProduct,
                                           geometryBuffer)
  { }

  /**
   * \brief Assemble the DPG system for a given rhs function.
   *
   * Given a LinearForm \p rhsLinearForm,
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
   * \tparam spaceIndex  the index of the solution space on which we apply
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
   * \brief The same as applyDirichletBoundary but it only
   *        applies the boundary values to the matrix
   *
   * \param[in,out] matrix      the matrix of the DPG system
   * \param[in] dirichletNodes  true marks the dofs in the Dirichlet boundary
   * \tparam spaceIndex  the index of the solution space on which we apply
   *                     the boundary data
   */
  template <size_t spaceIndex>
  void applyDirichletBoundaryToMatrix(
                              BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                              const std::vector<bool>& dirichletNodes);

  /**
   * \brief The same as applyDirichletBoundary but it only
   *        applies the boundary values to the rhs
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
   * \brief Apply Dirichlet boundary values to a solution space
   *
   * \param[in,out] matrix      the matrix of the DPG system
   * \param[in,out] rhs         the rhs vector of the DPG system
   * \param[in] dirichletNodes  true marks the dofs in the Dirichlet boundary
   * \param[in] value           the Dirichlet boundary value
   * \tparam spaceIndex  the index of the solution space on which we apply
   *                     the boundary data
   * \tparam ValueType   for functions for \p value
   */
  template <size_t spaceIndex, class ValueType>
  void applyNonzeroDirichletBoundary(
                              BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                              BlockVector<FieldVector<double,1> >& rhs,
                              const std::vector<bool>& dirichletNodes,
                              const ValueType& value);


  template <size_t spaceIndex, int dim>
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
  template <size_t spaceIndex, int dim>
  void defineCharacteristicFaces(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                                 BlockVector<FieldVector<double,1> >& rhs,
                                 const FieldVector<double,dim>& beta,
                                 double delta = 10e-10);

  template <size_t spaceIndex, class MinInnerProduct, int dim>
  void applyMinimization(
                      BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                      MinInnerProduct& minInnerProduct,
                      FieldVector<double, dim> beta,
                      double delta = 10e-10,
                      double epsilon = 0);

  /**
   * \brief Does exactly what it says on the tin.
   */
  TestSearchSpacesPtr getTestSearchSpaces() const
  { return testSearchSpaces_; }

  /**
   * \brief Does exactly what it says on the tin.
   */
  SolutionSpacesPtr getSolutionSpaces() const
  { return solutionSpaces_; }

private:
  template<size_t spaceIndex, int dim,
    typename std::enable_if<models<Functions::Concept::GlobalBasis<typename
                        std::tuple_element_t<spaceIndex,
                                             SolutionSpaces>::GridView>,
                  std::tuple_element_t<spaceIndex,
                                       SolutionSpaces>>()>::type* = nullptr>
  void defineCharacteristicFaces_impl(
      BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
      BlockVector<FieldVector<double,1> >& rhs,
      const FieldVector<double,dim>& beta,
      double delta);

  template<size_t spaceIndex, int dim,
    typename std::enable_if<models<Functions::Concept::ConstrainedGlobalBasis<
                      typename std::tuple_element_t<spaceIndex,
                                                    SolutionSpaces>::GridView>,
                  std::tuple_element_t<spaceIndex,
                                       SolutionSpaces>>()>::type* = nullptr>
  void defineCharacteristicFaces_impl(
      BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
      BlockVector<FieldVector<double,1> >& rhs,
      const FieldVector<double,dim>& beta,
      double delta);

  TestSearchSpacesPtr         testSearchSpaces_;
  SolutionSpacesPtr           solutionSpaces_;
  BilinearForm&               bilinearForm_;
  TestspaceCoefficientMatrix  testspaceCoefficientMatrix_;
};


/**
 * \brief Creates a SystemAssembler for a DPG formulation,
 *        deducing the target type from the types of arguments.
 *
 * \param bilinearForm   the bilinear form describing the DPG system
 * \param innerProduct   the inner product of the test search space
 */
template<class BilinearForm,
         class InnerProduct>
auto make_DPGSystemAssembler(BilinearForm&   bilinearForm,
                             InnerProduct&   innerProduct)
    -> DPGSystemAssembler<BilinearForm, InnerProduct, detail::Unbuffered>
{
  return DPGSystemAssembler<BilinearForm, InnerProduct, detail::Unbuffered>
                      (bilinearForm,
                       innerProduct);
}

/**
 * \brief Creates a SystemAssembler for a DPG formulation,
 *        deducing the target type from the types of arguments.
 *
 * \param bilinearForm   the bilinear form describing the DPG system
 * \param innerProduct   the inner product of the test search space
 * \param geometryBuffer
 */
template<class BilinearForm,
         class InnerProduct,
         class GeometryBuffer>
auto make_DPGSystemAssembler(BilinearForm&   bilinearForm,
                             InnerProduct&   innerProduct,
                             GeometryBuffer& geometryBuffer
                            )
    -> DPGSystemAssembler<BilinearForm, InnerProduct, detail::Buffered>
{
  return DPGSystemAssembler<BilinearForm, InnerProduct, detail::Buffered>
                      (bilinearForm,
                       innerProduct,
                       geometryBuffer);
}


template<class BilinearForm, class InnerProduct, class BufferPolicy>
template <class LinearForm>
void DPGSystemAssembler<BilinearForm, InnerProduct, BufferPolicy>::
assembleSystem(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
               BlockVector<FieldVector<double,1> >& rhs,
               LinearForm& rhsLinearForm)
{
  using namespace Dune::detail;

  const auto gridView = std::get<0>(*testSearchSpaces_).gridView();

  /* set up global offsets */
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

  const size_t globalTotalSolutionSize =
      computeOffsets(globalSolutionSpaceOffsets, *solutionSpaces_);

  // Views on the FE bases on a single element
  auto testLocalViews = getLocalViews(*testSearchSpaces_);

  auto solutionLocalViews = getLocalViews(*solutionSpaces_);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
  auto solutionLocalIndexSets = getLocalIndexSets(*solutionSpaces_);
#endif

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // TODO: Might be too large??
  MatrixIndexSet occupationPattern;
  occupationPattern.resize(globalTotalSolutionSize, globalTotalSolutionSize);

  namespace hana = boost::hana;
  auto indexRange = hana::make_range(hana::int_c<0>,
            hana::int_c<std::tuple_size<SolutionSpaces>::value>);
  auto indices = hana::cartesian_product(hana::make_tuple(indexRange,
                                                          indexRange));
  using Indices = decltype(indices);

  for(const auto& e : elements(gridView))
  {
    bindLocalViews(solutionLocalViews, e);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    bindLocalIndexSets(solutionLocalIndexSets, solutionLocalViews);
#endif

    detail::getOccupationPattern<Indices, false>
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
                        (solutionLocalViews,
                         solutionLocalViews,
                         globalSolutionSpaceOffsets,
                         globalSolutionSpaceOffsets,
                         occupationPattern);
#else
                        (solutionLocalIndexSets,
                         solutionLocalIndexSets,
                         globalSolutionSpaceOffsets,
                         globalSolutionSpaceOffsets,
                         occupationPattern);
#endif
  }

  occupationPattern.exportIdx(matrix);

  // set rhs to correct length -- the total number of basis vectors in the basis
  rhs.resize(globalTotalSolutionSize);

  // Set all entries to zero
  matrix = 0;
  rhs = 0;

  for(const auto& e : elements(gridView)) {

    bindLocalViews(solutionLocalViews, e);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    bindLocalIndexSets(solutionLocalIndexSets, solutionLocalViews);
#endif

    size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
    computeOffsets(localSolutionSpaceOffsets, solutionLocalViews);

    // compute the coefficient matrix C for the optimal test space
    testspaceCoefficientMatrix_.bind(e);
    const Matrix<FieldMatrix<double,1,1> >& coefficientMatrix
        = testspaceCoefficientMatrix_.coefficientMatrix();

    // Now get the local contribution to the right-hand side vector

    bindLocalViews(testLocalViews, e);

    // compute the local right-hand side vector F for the enriched test space
    BlockVector<FieldVector<double,1> > localEnrichedRhs;

    rhsLinearForm.bind(testLocalViews);
    rhsLinearForm.getLocalVector(localEnrichedRhs);

    // compute the local right-hand side vector C^T*F for the optimal test space
    BlockVector<FieldVector<double,1> > localRhs;
    localRhs.resize(coefficientMatrix.M());
    for (unsigned int i=0; i<coefficientMatrix.M(); i++)
      {
        localRhs[i]=0;
        for (unsigned int k=0; k<coefficientMatrix.N(); k++)
        {
          localRhs[i]+=(localEnrichedRhs[k]*coefficientMatrix[k][i]);
        }
      }

    // compute the local stiffness matrix
    const Matrix<FieldMatrix<double,1,1> >& elementMatrix
        = testspaceCoefficientMatrix_.systemMatrix();

    // Add element stiffness matrix onto the global stiffness matrix
    copyLocalToGlobalMatrix<Indices>(
        elementMatrix,
        matrix,
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
        solutionLocalViews,
#else
        solutionLocalIndexSets,
#endif
        localSolutionSpaceOffsets,
        globalSolutionSpaceOffsets,
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
        solutionLocalViews,
#else
        solutionLocalIndexSets,
#endif
        localSolutionSpaceOffsets,
        globalSolutionSpaceOffsets);

    // TODO: shouldn't LFIndices be taken from rhsLinearForm.getTerms()?
    namespace hana = boost::hana;
    using BilinearTerms = std::decay_t<decltype(bilinearForm_.getTerms())>;
    auto lfIndices = hana::to<hana::set_tag>(
        hana::transform(hana::to<hana::tuple_tag>(
            hana::make_range(hana::int_c<0>,
                hana::int_c<std::tuple_size<BilinearTerms>::value>)),
          [](auto i) -> auto {
            using Term = std::tuple_element_t<i.value, BilinearTerms>;
            return hana::type_c<typename Term::RhsIndex>;
          }));
    using LFIndices = decltype(hana::to<hana::tuple_tag>(lfIndices));

    // Add local right-hand side onto the global right-hand side
    copyLocalToGlobalVector<LFIndices>(
        localRhs,
        rhs,
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
        solutionLocalViews,
#else
        solutionLocalIndexSets,
#endif
        localSolutionSpaceOffsets,
        globalSolutionSpaceOffsets);
  }
}


template<class BilinearForm, class InnerProduct, class BufferPolicy>
void DPGSystemAssembler<BilinearForm, InnerProduct, BufferPolicy>::
assembleMatrix(BCRSMatrix<FieldMatrix<double,1,1> >& matrix)
{
  using namespace Dune::detail;

  const auto gridView = std::get<0>(*testSearchSpaces_).gridView();

  /* set up global offsets */
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

  const size_t globalTotalSolutionSize =
      computeOffsets(globalSolutionSpaceOffsets, *solutionSpaces_);

  // Views on the FE bases on a single element
  auto solutionLocalViews = getLocalViews(*solutionSpaces_);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
  auto solutionLocalIndexSets = getLocalIndexSets(*solutionSpaces_);
#endif

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // TODO: Might be too large??
  MatrixIndexSet occupationPattern;
  occupationPattern.resize(globalTotalSolutionSize, globalTotalSolutionSize);

  namespace hana = boost::hana;
  auto indexRange = hana::make_range(hana::int_c<0>,
            hana::int_c<std::tuple_size<SolutionSpaces>::value>);
  auto indices = hana::cartesian_product(hana::make_tuple(indexRange,
                                                          indexRange));
  using Indices = decltype(indices);

  for(const auto& e : elements(gridView))
  {
    bindLocalViews(solutionLocalViews, e);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    bindLocalIndexSets(solutionLocalIndexSets, solutionLocalViews);
#endif

    detail::getOccupationPattern<Indices, false>(
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
                         solutionLocalViews,
                         solutionLocalViews,
#else
                         solutionLocalIndexSets,
                         solutionLocalIndexSets,
#endif
                         globalSolutionSpaceOffsets,
                         globalSolutionSpaceOffsets,
                         occupationPattern);
  }

  occupationPattern.exportIdx(matrix);

  // Set all entries to zero
  matrix = 0;

  for(const auto& e : elements(gridView)) {

    bindLocalViews(solutionLocalViews, e);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    bindLocalIndexSets(solutionLocalIndexSets, solutionLocalViews);
#endif

    testspaceCoefficientMatrix_.bind(e);
    const Matrix<FieldMatrix<double,1,1> >& elementMatrix
        = testspaceCoefficientMatrix_.systemMatrix();

    size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
    computeOffsets(localSolutionSpaceOffsets, solutionLocalViews);

    // Add element stiffness matrix onto the global stiffness matrix
    /* copy every local submatrix indexed by a pair of indices from
     * Indices exactly once. */
    copyLocalToGlobalMatrix<Indices>(
        elementMatrix,
        matrix,
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
        solutionLocalViews,
#else
        solutionLocalIndexSets,
#endif
        localSolutionSpaceOffsets,
        globalSolutionSpaceOffsets,
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
        solutionLocalViews,
#else
        solutionLocalIndexSets,
#endif
        localSolutionSpaceOffsets,
        globalSolutionSpaceOffsets);
  }
}


template<class BilinearForm, class InnerProduct, class BufferPolicy>
template <class LinearForm>
void DPGSystemAssembler<BilinearForm, InnerProduct, BufferPolicy>::
assembleRhs(BlockVector<FieldVector<double,1> >& rhs,
            LinearForm& rhsLinearForm)
{
  using namespace Dune::detail;

  const auto gridView = std::get<0>(*testSearchSpaces_).gridView();

  /* set up global offsets */
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

  const size_t globalTotalSolutionSize =
      computeOffsets(globalSolutionSpaceOffsets, *solutionSpaces_);

  // set rhs to correct length -- the total number of basis vectors in the basis
  rhs.resize(globalTotalSolutionSize);

  // Set all entries to zero
  rhs = 0;

  // Views on the FE bases on a single element
  auto testLocalViews = getLocalViews(*testSearchSpaces_);

  auto solutionLocalViews = getLocalViews(*solutionSpaces_);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
  auto solutionLocalIndexSets = getLocalIndexSets(*solutionSpaces_);
#endif

  for(const auto& e : elements(gridView)) {

    bindLocalViews(solutionLocalViews, e);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    bindLocalIndexSets(solutionLocalIndexSets, solutionLocalViews);
#endif

    size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
    computeOffsets(localSolutionSpaceOffsets, solutionLocalViews);

    // compute the coefficient matrix C for the optimal test space
    testspaceCoefficientMatrix_.bind(e);
    const Matrix<FieldMatrix<double,1,1> >& coefficientMatrix
        = testspaceCoefficientMatrix_.coefficientMatrix();

    // Now get the local contribution to the right-hand side vector

    bindLocalViews(testLocalViews, e);

    // compute the local right-hand side vector F for the enriched test space
    BlockVector<FieldVector<double,1> > localEnrichedRhs;

    rhsLinearForm.bind(testLocalViews);
    rhsLinearForm.getLocalVector(localEnrichedRhs);

    // compute the local right-hand side vector C^T*F for the optimal test space
    BlockVector<FieldVector<double,1> > localRhs;
    localRhs.resize(coefficientMatrix.M());
    for (unsigned int i=0; i<coefficientMatrix.M(); i++)
      {
        localRhs[i]=0;
        for (unsigned int k=0; k<coefficientMatrix.N(); k++)
        {
          localRhs[i]+=(localEnrichedRhs[k]*coefficientMatrix[k][i]);
        }
      }

    /* copy every local subvector indexed by an index from
     * lfIndices exactly once. */
    namespace hana = boost::hana;
    // TODO: shouldn't LFIndices be taken from rhsLinearForm.getTerms()?
    using BilinearTerms = std::decay_t<decltype(bilinearForm_.getTerms())>;
    auto lfIndices = hana::to<hana::set_tag>(
        hana::transform(hana::make_range(hana::int_c<0>,
              hana::int_c<std::tuple_size<BilinearTerms>::value>),
          [](auto i) {
            using Term = std::tuple_element_t<i.value, BilinearTerms>;
            return hana::type_c<typename Term::RhsIndex>;
          }));
    using LFIndices = decltype(lfIndices);

    // Add local right-hand side onto the global right-hand side
    copyLocalToGlobalVector<LFIndices>(
        localRhs,
        rhs,
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
        solutionLocalViews,
#else
        solutionLocalIndexSets,
#endif
        localSolutionSpaceOffsets,
        globalSolutionSpaceOffsets);

  }
}


template<class BilinearForm, class InnerProduct, class BufferPolicy>
template <size_t spaceIndex, class ValueType>
void DPGSystemAssembler<BilinearForm, InnerProduct, BufferPolicy>::
applyDirichletBoundary
                      (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                       BlockVector<FieldVector<double,1> >& rhs,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& boundaryValue)
{
  applyDirichletBoundaryToMatrix<spaceIndex>
          (matrix, dirichletNodes);
  applyDirichletBoundaryToRhs<spaceIndex>
          (rhs,    dirichletNodes, boundaryValue);
}


template<class BilinearForm, class InnerProduct, class BufferPolicy>
template <size_t spaceIndex>
void DPGSystemAssembler<BilinearForm, InnerProduct, BufferPolicy>::
applyDirichletBoundaryToMatrix
                      (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                       const std::vector<bool>& dirichletNodes)
{
  const size_t spaceSize =
        std::get<spaceIndex>(*solutionSpaces_).size();

  const size_t globalOffset =
        detail::computeOffset<spaceIndex>(*solutionSpaces_);

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


template<class BilinearForm, class InnerProduct, class BufferPolicy>
template <size_t spaceIndex, class ValueType>
void DPGSystemAssembler<BilinearForm, InnerProduct, BufferPolicy>::
applyDirichletBoundaryToRhs
                      (BlockVector<FieldVector<double,1> >& rhs,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& boundaryValue)
{
  const size_t spaceSize =
        std::get<spaceIndex>(*solutionSpaces_).size();

  const size_t globalOffset =
        detail::computeOffset<spaceIndex>(*solutionSpaces_);

  // Set Dirichlet values
  for (size_t i=0; i<spaceSize; i++)
  {
    if (dirichletNodes[i])
    {
      /* TODO: Needs adaptation when value is a function. */
      rhs[globalOffset+i] = detail::evaluateBoundary(boundaryValue,i);
    }
  }

}


template<class BilinearForm, class InnerProduct, class BufferPolicy>
template <size_t spaceIndex, class ValueType>
void DPGSystemAssembler<BilinearForm, InnerProduct, BufferPolicy>::
applyNonzeroDirichletBoundary(
                              BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                              BlockVector<FieldVector<double,1> >& rhs,
                              const std::vector<bool>& dirichletNodes,
                              const ValueType& value)
{
  const size_t spaceSize =
        std::get<spaceIndex>(*solutionSpaces_).size();

  const size_t globalOffset =
        detail::computeOffset<spaceIndex>(*solutionSpaces_);

  ////////////////////////////////////////////
  //    Modify rhs vector
  ////////////////////////////////////////////

  // compute the coefficients for the Dirichlet nodes
  std::vector<double> dirichletValues;
  dirichletValues.resize(std::get<spaceIndex>(*solutionSpaces_).size());
  interpolate(std::get<spaceIndex>(*solutionSpaces_),
              Dune::TypeTree::hybridTreePath(),
              dirichletValues, value,
              dirichletNodes);

  ////////////////////////////////////////////
  //   Modify Dirichlet rows in matrix
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
          /* Zero out row and column to keep symmetry. Modify RHS accordingly*/
          *cIt = 0;
          rhs[cIt.index()]-=(matrix[cIt.index()][globalOffset+i]*dirichletValues[i]);
          matrix[cIt.index()][globalOffset+i]=0;
        }
      }
    }

  }

  // Set Dirichlet values in rhs
  for (size_t i=0; i<spaceSize; i++)
  {
    if (dirichletNodes[i])
    {
      rhs[globalOffset+i] = dirichletValues[i];
    }
  }
}





template<class BilinearForm, class InnerProduct, class BufferPolicy>
template <size_t spaceIndex, int dim>
void DPGSystemAssembler<BilinearForm, InnerProduct, BufferPolicy>::
applyWeakBoundaryCondition
                    (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                     const FieldVector<double, dim> beta,
                     const double mu)
{
  const auto gridView = std::get<0>(*solutionSpaces_).gridView();

  const size_t globalOffset =
        detail::computeOffset<spaceIndex>(*solutionSpaces_);

  auto localView     = std::get<spaceIndex>(*solutionSpaces_).localView();
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
  auto localIndexSet = std::get<spaceIndex>(*solutionSpaces_).localIndexSet();
#endif

  for(const auto& e : elements(gridView))
  {
    localView.bind(e);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    localIndexSet.bind(localView);
#endif

    const auto& localFiniteElement = localView.tree().finiteElement();

    const unsigned int quadratureOrder
        = 2*localFiniteElement.localBasis().order();

    const size_t n = localFiniteElement.size();

    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    // Set all matrix entries to zero
    elementMatrix.setSize(n,n);
    elementMatrix = 0;

    for (auto&& intersection : intersections(gridView, e))
    {
      if (intersection.boundary())
      { // the intersection is at the (physical) boundary of the domain
        const FieldVector<double,dim>& centerOuterNormal =
               intersection.centerUnitOuterNormal();

        if ((beta*centerOuterNormal) > -1e-10)
        { // everywhere except inflow boundary
          const QuadratureRule<double, dim-1>& quadFace =
                  QuadratureRules<double, dim-1>::rule(intersection.type(),
                                                       quadratureOrder);

          for (size_t pt=0; pt < quadFace.size(); pt++)
          {
            // position of the current quadrature point in the
            // reference element (face!)
            const FieldVector<double,dim-1>& quadFacePos
                = quadFace[pt].position();

            const double integrationWeight
                = intersection.geometry().integrationElement(quadFacePos)
                * quadFace[pt].weight();

            // position of the quadrature point within the element
            const FieldVector<double,dim> elementQuadPos
                = intersection.geometryInInside().global(quadFacePos);

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
                           * integrationWeight;
              }
            }

          }
        }
      }
    }
    addToGlobalMatrix(
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
        localView,
        localView,
#else
        localIndexSet,
        localIndexSet,
#endif
        [&elementMatrix](size_t i, size_t j) -> auto {
          return elementMatrix[i][j];
        },
        [&](auto gi, auto gj) -> auto& {
          return matrix[gi[0]+globalOffset][gj[0]+globalOffset];
        }
    );
  }
}


template<class BilinearForm, class InnerProduct, class BufferPolicy>
template<size_t spaceIndex, int dim,
  typename std::enable_if<models<Functions::Concept::GlobalBasis<typename
                       std::tuple_element_t<spaceIndex, typename
                          BilinearForm::SolutionSpaces>::GridView>,
                 std::tuple_element_t<spaceIndex, typename
                       BilinearForm::SolutionSpaces>>()>::type*>
void DPGSystemAssembler<BilinearForm, InnerProduct, BufferPolicy>::
defineCharacteristicFaces_impl(
    BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
    BlockVector<FieldVector<double,1> >& rhs,
    const FieldVector<double,dim>& beta,
    const double delta)
{
  const auto gridView = std::get<spaceIndex>(*solutionSpaces_).gridView();

  const size_t globalOffset =
        detail::computeOffset<spaceIndex>(*solutionSpaces_);

  auto solutionLocalView = std::get<spaceIndex>(*solutionSpaces_).localView();
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
  auto localIndexSet = std::get<spaceIndex>(*solutionSpaces_).localIndexSet();
#endif

  for(const auto& e : elements(gridView))
  {
    solutionLocalView.bind(e);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    localIndexSet.bind(solutionLocalView);
#endif

    const auto& localFiniteElement = solutionLocalView.tree().finiteElement();

    const size_t n = localFiniteElement.size();

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
      if(e.type().isQuadrilateral()) {
        endpoints[0] = std::make_pair(vertexDOFs[0], vertexDOFs[2]);
        endpoints[1] = std::make_pair(vertexDOFs[1], vertexDOFs[3]);
        endpoints[2] = std::make_pair(vertexDOFs[0], vertexDOFs[1]);
        endpoints[3] = std::make_pair(vertexDOFs[2], vertexDOFs[3]);
      } else if(e.type().isTriangle()) {
        endpoints[0] = std::make_pair(vertexDOFs[0], vertexDOFs[1]);
        endpoints[1] = std::make_pair(vertexDOFs[0], vertexDOFs[2]);
        endpoints[2] = std::make_pair(vertexDOFs[1], vertexDOFs[2]);
      } else {
        DUNE_THROW(Dune::NotImplemented,
                   "defineCharacteristicFaces not implemented for element type"
                   << e.type().id());
      }

      for (auto&& faceAndDOFs: characteristicDOFs)
      {
        size_t face;
        std::list<std::pair<size_t,size_t>> dofs;
        std::tie(face, dofs) = faceAndDOFs;
        size_t left, right;
        std::tie(left, right) = endpoints[face];
        for(auto&& dof: dofs)
        {
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
          const auto row = solutionLocalView.index(dof.first)[0];
#else
          const auto row = localIndexSet.index(dof.first)[0];
#endif
          auto col = row;
          const size_t k = dofs.size()+1;

          /* replace the row of dof on characteristic face
           * by an interpolation of the two endpoints of the
           * characteristic face. */
          matrix[row+globalOffset][col+globalOffset] = -1;
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
          col = solutionLocalView.index(left)[0];
#else
          col = localIndexSet.index(left)[0];
#endif
          matrix[row+globalOffset][col+globalOffset]
              = static_cast<double>(k-dof.second-1)/k;
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
          col = solutionLocalView.index(right)[0];
#else
          col = localIndexSet.index(right)[0];
#endif
          matrix[row+globalOffset][col+globalOffset]
              = static_cast<double>(dof.second+1)/k;

          rhs[row+globalOffset] = 0;
        }
      }
    }
  }
}


namespace detail {
  template<class ConstrainedLocalIndexSet>
  typename ConstrainedLocalIndexSet::MultiIndex
  getUnconstrainedIndex(const ConstrainedLocalIndexSet& localIndexSet,
      size_t localIndex)
  {
    using size_type
        = typename std::decay_t<ConstrainedLocalIndexSet>::size_type;
    auto globalIndex = localIndexSet.indicesLocalGlobal().begin();
    const size_type numConstraints = localIndexSet.constraintsSize();
    size_type i = 0;
    for(size_type c = 0; c < numConstraints; c++) {
      const size_type nextConstraint = i + localIndexSet.constraintOffset(c);
      if(localIndex < nextConstraint) {
        return globalIndex[localIndex - i];
      } else {
        assert(localIndex != nextConstraint);
        std::advance(globalIndex, localIndexSet.constraintOffset(c)
            + localIndexSet.constraintWeights(c).size());
        i = nextConstraint + 1;
      }
    }
    return globalIndex[localIndex - i];
  }
}


template<class BilinearForm, class InnerProduct, class BufferPolicy>
template<size_t spaceIndex, int dim,
  typename std::enable_if<models<Functions::Concept::ConstrainedGlobalBasis<
                    typename std::tuple_element_t<spaceIndex, typename
                        BilinearForm::SolutionSpaces>::GridView>,
                std::tuple_element_t<spaceIndex, typename
                      BilinearForm::SolutionSpaces>>()>::type*>
void DPGSystemAssembler<BilinearForm, InnerProduct, BufferPolicy>::
defineCharacteristicFaces_impl(
    BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
    BlockVector<FieldVector<double,1> >& rhs,
    const FieldVector<double,dim>& beta,
    const double delta)
{
  const auto gridView = std::get<spaceIndex>(*solutionSpaces_).gridView();

  const size_t globalOffset =
        detail::computeOffset<spaceIndex>(*solutionSpaces_);

  auto solutionLocalView = std::get<spaceIndex>(*solutionSpaces_).localView();
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
  auto localIndexSet = std::get<spaceIndex>(*solutionSpaces_).localIndexSet();
#endif

  for(const auto& e : elements(gridView))
  {
    solutionLocalView.bind(e);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    localIndexSet.bind(solutionLocalView);
#endif

    const auto& localFiniteElement = solutionLocalView.tree().finiteElement();

    const size_t n = localFiniteElement.size();

    std::vector<bool>  characteristicFaces(e.subEntities(1), false);
    bool characteristicFound = false;

    for (auto&& intersection : intersections(gridView, e))
    {
      if (intersection.conforming()) {
        const bool characteristic =
            fabs(beta * intersection.centerUnitOuterNormal()) < delta;
        characteristicFaces[intersection.indexInInside()] = characteristic;
        characteristicFound = characteristicFound || characteristic;
      }
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
      if(e.type().isQuadrilateral()) {
        endpoints[0] = std::make_pair(vertexDOFs[0], vertexDOFs[2]);
        endpoints[1] = std::make_pair(vertexDOFs[1], vertexDOFs[3]);
        endpoints[2] = std::make_pair(vertexDOFs[0], vertexDOFs[1]);
        endpoints[3] = std::make_pair(vertexDOFs[2], vertexDOFs[3]);
      } else if(e.type().isTriangle()) {
        endpoints[0] = std::make_pair(vertexDOFs[0], vertexDOFs[1]);
        endpoints[1] = std::make_pair(vertexDOFs[0], vertexDOFs[2]);
        endpoints[2] = std::make_pair(vertexDOFs[1], vertexDOFs[2]);
      } else {
        DUNE_THROW(Dune::NotImplemented,
                   "defineCharacteristicFaces not implemented for element type"
                   << e.type().id());
      }

      for (auto&& faceAndDOFs: characteristicDOFs)
      {
        size_t face;
        std::list<std::pair<size_t,size_t>> dofs;
        std::tie(face, dofs) = faceAndDOFs;
        size_t left, right;
        std::tie(left, right) = endpoints[face];
        for(auto&& dof: dofs)
        {
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
          const auto row = detail::getUnconstrainedIndex(solutionLocalView,
                                                         dof.first)[0];
#else
          const auto row = detail::getUnconstrainedIndex(localIndexSet,
                                                         dof.first)[0];
#endif
          auto col = row;
          const size_t k = dofs.size()+1;

          /* replace the row of dof on characteristic face
           * by an interpolation of the two endpoints of the
           * characteristic face. */
          matrix[row+globalOffset][col+globalOffset] = -1;
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
          col = detail::getUnconstrainedIndex(solutionLocalView, left)[0];
#else
          col = detail::getUnconstrainedIndex(localIndexSet, left)[0];
#endif
          matrix[row+globalOffset][col+globalOffset]
              = static_cast<double>(k-dof.second-1)/k;
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
          col = detail::getUnconstrainedIndex(solutionLocalView, right)[0];
#else
          col = detail::getUnconstrainedIndex(localIndexSet, right)[0];
#endif
          matrix[row+globalOffset][col+globalOffset]
              = static_cast<double>(dof.second+1)/k;

          rhs[row+globalOffset] = 0;
        }
      }
    }
  }
}


template<class BilinearForm, class InnerProduct, class BufferPolicy>
template<size_t spaceIndex, int dim>
void DPGSystemAssembler<BilinearForm, InnerProduct, BufferPolicy>::
defineCharacteristicFaces(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                          BlockVector<FieldVector<double,1> >& rhs,
                          const FieldVector<double,dim>& beta,
                          const double delta)
{
  // make sure that the test search spaces are unrefined
  {
    namespace hana = boost::hana;
    auto spacesRefined =
        hana::transform(hana::to<hana::tuple_tag>(
            hana::make_range(hana::int_c<0>,
                hana::int_c<std::tuple_size<TestSearchSpaces>::value>)),
          [](auto i) -> auto {
            using Space = std::tuple_element_t<i.value, TestSearchSpaces>;
            using Refined = typename is_RefinedFiniteElement<Space>::type;
            return Refined{};
          });
    static_assert(hana::none(spacesRefined),
        "defineCharacteristicFaces only works for unrefined test spaces.");
  }
  defineCharacteristicFaces_impl<spaceIndex, dim>(
      matrix,
      rhs,
      beta,
      delta);
}


template<class BilinearForm, class InnerProduct, class BufferPolicy>
template <size_t spaceIndex, class MinInnerProduct, int dim>
void DPGSystemAssembler<BilinearForm, InnerProduct, BufferPolicy>::
applyMinimization
            (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
             MinInnerProduct& minInnerProduct,
             const FieldVector<double, dim> beta,
             const double delta,
             const double epsilon)
{
  using namespace Dune::detail;

  const auto gridView = std::get<spaceIndex>(*solutionSpaces_).gridView();

  //const size_t globalOffset = computeOffset<spaceIndex>(*solutionSpaces_);
  const size_t globalOffset = 0;

  size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

  // get local view for solution space
  // (necessary if we want to use inner product) // TODO inefficient (why?)
  auto solutionLocalViews = getLocalViews(*solutionSpaces_);

#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
  auto localIndexSet = std::get<spaceIndex>(*solutionSpaces_).localIndexSet();
  using LocalIndexSet = decltype(localIndexSet);
#else
  using LocalView = std::tuple_element_t<spaceIndex,
                                         decltype(solutionLocalViews)>;
#endif

  const bool epsilonSmallerDelta = epsilon<delta;

  for(const auto& e : elements(gridView))
  {
    bindLocalViews(solutionLocalViews, e);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    localIndexSet.bind(std::get<spaceIndex>(solutionLocalViews));
#endif

    /* set up local offsets */
    computeOffsets(localSolutionSpaceOffsets, solutionLocalViews);

    const auto& localFiniteElement
        = std::get<spaceIndex>(solutionLocalViews).tree().finiteElement();

    const size_t n = localFiniteElement.size();

    Matrix<FieldMatrix<double,1,1> > elementMatrix;

    minInnerProduct.bind(solutionLocalViews);
    minInnerProduct.getLocalMatrix(elementMatrix);

    std::vector<bool> relevantFaces(e.subEntities(1), false);

    for (auto&& intersection : intersections(gridView, e))
    {
      const FieldVector<double,dim>& centerOuterNormal =
             intersection.centerUnitOuterNormal();
      //set relevant faces for almost characteristic faces
      relevantFaces[intersection.indexInInside()]
          = (std::abs(beta*centerOuterNormal) < delta);
    }

    std::vector<bool> relevantDOFs(n, false);

    for (unsigned int i=0; i<n; i++)
    {
      if (localFiniteElement.localCoefficients().localKey(i).codim()==0)
      { // interior DOFs
        relevantDOFs[i] = true;
      }
      else if (localFiniteElement.localCoefficients().localKey(i).codim()==1
               and epsilonSmallerDelta)
      { // edge DOFs
        relevantDOFs[i] = relevantFaces[localFiniteElement.localCoefficients()
                                          .localKey(i).subEntity()];
      }
      // Vertex DOFs are never relevant because the corresponding
      // basis functions have support on at least two edges which can
      // never be both (almost) characteristic.
    }

#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
    using MultiIndex = typename std::decay_t<LocalView>::MultiIndex;
#else
    using MultiIndex = typename std::decay_t<LocalIndexSet>::MultiIndex;
#endif
    iterateOverLocalIndices(
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
      std::get<spaceIndex>(solutionLocalViews),
#else
      localIndexSet,
#endif
      [&](size_t i, MultiIndex gi)
      {
        if (relevantDOFs[i])
        {
          auto row = gi[0];
          iterateOverLocalIndices(
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
            std::get<spaceIndex>(solutionLocalViews),
#else
            localIndexSet,
#endif
            [&](size_t j, MultiIndex gj)
            {
              auto col = gj[0];
              matrix[row+globalOffset][col+globalOffset]
                  += elementMatrix[localSolutionSpaceOffsets[spaceIndex]+i]
                                  [localSolutionSpaceOffsets[spaceIndex]+j];
            },
            [](size_t j) {},
            [&](size_t j, MultiIndex gj, double wj)
            {
              auto col = gj[0];
              matrix[row+globalOffset][col+globalOffset]
                += wj * elementMatrix[localSolutionSpaceOffsets[spaceIndex]+i]
                                     [localSolutionSpaceOffsets[spaceIndex]+j];
            }
          );
        }
      },
      [](size_t i) {},
      [&](size_t i, MultiIndex gi, double wi)
      {
        if (relevantDOFs[i])
        {
          auto row = gi[0];
          iterateOverLocalIndices(
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
            std::get<spaceIndex>(solutionLocalViews),
#else
            localIndexSet,
#endif
            [&](size_t j, MultiIndex gj)
            {
              auto col = gj[0];
              matrix[row+globalOffset][col+globalOffset]
                  += wi
                     * elementMatrix[localSolutionSpaceOffsets[spaceIndex]+i]
                                    [localSolutionSpaceOffsets[spaceIndex]+j];
            },
            [](size_t j) {},
            [&](size_t j, MultiIndex gj, double wj)
            {
              auto col = gj[0];
              matrix[row+globalOffset][col+globalOffset]
                += wi * wj
                   * elementMatrix[localSolutionSpaceOffsets[spaceIndex]+i]
                                  [localSolutionSpaceOffsets[spaceIndex]+j];
            }
          );
        }
      }
    );
  }
}

} // end namespace Dune

#endif // DUNE_DPG_SYSTEM_ASSEMBLER_DPG_HH
