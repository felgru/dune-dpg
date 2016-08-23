// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_SYSTEM_ASSEMBLER_DPG_HH
#define DUNE_DPG_SYSTEM_ASSEMBLER_DPG_HH

#include <functional>
#include <list>
#include <map>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

#include <boost/mpl/empty_sequence.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/joint_view.hpp>
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
#include "testspace_coefficient_matrix.hh"


namespace Dune {

/**
 * \brief This constructs the system matrix and righthandside vector of a DPG system.
 *
 * \tparam InnProduct        inner product
 * \tparam BilinForm         bilinear form
 */
template<class InnProduct, class BilinForm>
class DPGSystemAssembler
{
public:
  using InnerProduct = InnProduct;
  using BilinearForm = BilinForm;
  using TestspaceCoefficientMatrix = Dune::TestspaceCoefficientMatrix<BilinearForm, InnerProduct>;
  using TestSearchSpaces = typename BilinearForm::TestSpaces;
  using SolutionSpaces = typename BilinearForm::SolutionSpaces;


  //! tuple type for the local views of the test spaces
  typedef typename boost::fusion::result_of::as_vector<
      typename boost::fusion::result_of::
      transform<TestSearchSpaces, detail::getLocalView>::type>::type TestLocalViews;
  //! tuple type for the local views of the solution spaces
  typedef typename boost::fusion::result_of::as_vector<
      typename boost::fusion::result_of::
      transform<SolutionSpaces, detail::getLocalView>::type
      >::type SolutionLocalViews;

  DPGSystemAssembler () = delete;
  /**
   * \brief constructor for DPGSystemAssembler
   *
   * \note For your convenience, use make_DPGSystemAssembler()
   *       instead.
   */
  constexpr DPGSystemAssembler (InnerProduct&      innerProduct,
                                BilinearForm&      bilinearForm)
             : testSearchSpaces_(bilinearForm.getTestSpaces()),
               solutionSpaces_(bilinearForm.getSolutionSpaces()),
               bilinearForm_(bilinearForm),
               testspaceCoefficientMatrix_(bilinearForm,innerProduct)
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
   * \param[in] value           the Dirichlet boundary value
   * \tparam spaceIndex  the index of the solution space on which we apply
   *                     the boundary data
   * \tparam ValueType   we take either constants or functions for \p value
   */
  template <size_t spaceIndex, class ValueType>
  void applyDirichletBoundaryToMatrix(
                              BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                              const std::vector<bool>& dirichletNodes,
                              const ValueType& value);

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
  const TestSearchSpaces& getTestSearchSpaces() const
  { return testSearchSpaces_; }

  /**
   * \brief Does exactly what it says on the tin.
   */
  const SolutionSpaces& getSolutionSpaces() const
  { return solutionSpaces_; }

private:
  TestSearchSpaces              testSearchSpaces_;
  SolutionSpaces                solutionSpaces_;
  BilinearForm&                 bilinearForm_;
  TestspaceCoefficientMatrix    testspaceCoefficientMatrix_;
};


/**
 * \brief Creates a SystemAssembler for a DPG formulation,
 *        deducing the target type from the types of arguments.
 *
 * \param testSearchSpaces     a tuple of test spaces
 * \param solutionSpaces a tuple of solution spaces
 * \param bilinearForm   the bilinear form describing the DPG system
 */
template<class InnerProduct,
         class BilinearForm>
auto make_DPGSystemAssembler(InnerProduct&   innerProduct,
                             BilinearForm&   bilinearForm)
    -> DPGSystemAssembler<InnerProduct, BilinearForm>
{
  // set the inner product of the system assembler to nullptr as it is
  // not used in the DPG formulation (we use the inner product of the
  // optimal test space her).
  return DPGSystemAssembler<InnerProduct, BilinearForm>
                      (innerProduct,
                       bilinearForm);
}


template<class InnerProduct, class BilinearForm>
template <class LinearForm>
void DPGSystemAssembler<InnerProduct, BilinearForm>::
assembleSystem(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
               BlockVector<FieldVector<double,1> >& rhs,
               LinearForm& rhsLinearForm)
{
  assembleMatrix(matrix);
  assembleRhs   (rhs, rhsLinearForm);
}


template<class InnerProduct, class BilinearForm>
void DPGSystemAssembler<InnerProduct, BilinearForm>::
assembleMatrix(BCRSMatrix<FieldMatrix<double,1,1> >& matrix)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  typedef typename std::tuple_element<0,TestSearchSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSearchSpaces_).gridView();

  /* set up global offsets */
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

  size_t globalTotalSolutionSize =
      fold(zip(globalSolutionSpaceOffsets, solutionSpaces_),
           (size_t)0, globalOffsetHelper());

  globalTotalSolutionSize -= globalSolutionSpaceOffsets[0];

  const auto n = globalTotalSolutionSize;


  // Views on the FE bases on a single element
  auto solutionLocalViews = as_vector(transform(solutionSpaces_,
                                                getLocalView()));
  auto solutionLocalIndexSets = as_vector(transform(solutionSpaces_,
                                                    getLocalIndexSet()));

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // TODO: Might be too large??
  MatrixIndexSet occupationPattern;
  occupationPattern.resize(n, n);

  typedef
      typename result_of::as_vector<typename boost::mpl::range_c<
                            size_t,0,result_of::size<SolutionSpaces>
                         ::type::value>::type
            >::type IndexRange;
  typedef
      typename boost::mpl::fold<
            typename boost::mpl::transform<
                IndexRange
              , mpl::prefixSequenceWith<IndexRange, boost::mpl::_1>
              >::type
          , boost::mpl::empty_sequence
          , boost::mpl::joint_view<boost::mpl::_1, boost::mpl::_2>
          >::type Indices;

  for(const auto& e : elements(gridView))
  {
    for_each(solutionLocalViews, applyBind<decltype(e)>(e));

    for_each(zip(solutionLocalIndexSets, solutionLocalViews),
             make_fused_procedure(bindLocalIndexSet()));

    auto gOPH = getOccupationPatternHelper<decltype(solutionLocalViews),
                                           decltype(solutionLocalViews),
                                           decltype(solutionLocalIndexSets),
                                           decltype(solutionLocalIndexSets),
                                           false>
                        (solutionLocalViews,
                         solutionLocalViews,
                         solutionLocalIndexSets,
                         solutionLocalIndexSets,
                         globalSolutionSpaceOffsets,
                         globalSolutionSpaceOffsets,
                         occupationPattern);
    for_each(Indices{},
        std::ref(gOPH));
  }

  occupationPattern.exportIdx(matrix);

  // Set all entries to zero
  matrix = 0;

  for(const auto& e : elements(gridView)) {

    testspaceCoefficientMatrix_.bind(e);
    Matrix<FieldMatrix<double,1,1> > elementMatrix = testspaceCoefficientMatrix_.localMatrix();

    for_each(solutionLocalViews, applyBind<decltype(e)>(e));
    for_each(zip(solutionLocalIndexSets, solutionLocalViews),
             make_fused_procedure(bindLocalIndexSet()));

    // Add element stiffness matrix onto the global stiffness matrix
    auto cp = fused_procedure<localToGlobalCopier<decltype(elementMatrix),
                   typename std::remove_reference<decltype(matrix)>::type> >
                (localToGlobalCopier<decltype(elementMatrix),
                   typename std::remove_reference<decltype(matrix)>::type>
                                        (elementMatrix, matrix));

    /* copy every local submatrix indexed by a pair of indices from
     * Indices exactly once. */    auto solutionZip = zip(solutionLocalViews,
                           solutionLocalIndexSets,
                           bilinearForm_.getLocalSolutionSpaceOffsets(),
                           globalSolutionSpaceOffsets);

    using SolutionZip = decltype(solutionZip);

    for_each(Indices{},
             localToGlobalCopyHelper<SolutionZip,
                                     SolutionZip,
                                     std::decay_t<decltype(cp)>>
                                    (solutionZip, solutionZip, cp));
  }
}


template<class InnerProduct, class BilinearForm>
template <class LinearForm>
void DPGSystemAssembler<InnerProduct, BilinearForm>::
assembleRhs(BlockVector<FieldVector<double,1> >& rhs,
            LinearForm& rhsLinearForm)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  typedef typename std::tuple_element<0,TestSearchSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSearchSpaces_).gridView();

  /* set up global offsets */
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

  size_t globalTotalSolutionSize =
      fold(zip(globalSolutionSpaceOffsets, solutionSpaces_),
           (size_t)0, globalOffsetHelper());

  // set rhs to correct length -- the total number of basis vectors in the basis
  rhs.resize(globalTotalSolutionSize);

  // Set all entries to zero
  rhs = 0;

  // Views on the FE bases on a single element
  auto testLocalViews     = as_vector(transform(testSearchSpaces_, getLocalView()));

  auto solutionLocalViews     = as_vector(transform(solutionSpaces_, getLocalView()));

  auto solutionLocalIndexSets = as_vector(transform(solutionSpaces_,
                                                getLocalIndexSet()));

  for(const auto& e : elements(gridView)) {

    for_each(solutionLocalViews, applyBind<decltype(e)>(e));
    for_each(zip(solutionLocalIndexSets, solutionLocalViews),
             make_fused_procedure(bindLocalIndexSet()));

    size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
    fold(zip(localSolutionSpaceOffsets, solutionLocalViews),
               (size_t)0, offsetHelper());

    for_each(testLocalViews, applyBind<decltype(e)>(e));

    // Now get the local contribution to the right-hand side vector

    // compute the local right-hand side vector F for the enriched test space
    BlockVector<FieldVector<double,1> > localEnrichedRhs;

    rhsLinearForm.bind(testLocalViews);
    rhsLinearForm.getLocalVector(localEnrichedRhs);

    // compute the coefficient matrix C for the optimal test space
    testspaceCoefficientMatrix_.bind(e);
    Matrix<FieldMatrix<double,1,1> > coefficientMatrix = testspaceCoefficientMatrix_.coefficientMatrix();

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

    auto cp = fused_procedure<localToGlobalRHSCopier<decltype(localRhs),
                   typename std::remove_reference<decltype(rhs)>::type> >
                (localToGlobalRHSCopier<decltype(localRhs),
                   typename std::remove_reference<decltype(rhs)>::type>
                                        (localRhs, rhs));

    /* copy every local subvector indexed by an index from
     * lfIndices exactly once. */
    auto solutionZip = zip(solutionLocalViews,
                       solutionLocalIndexSets,
                       localSolutionSpaceOffsets,
                       globalSolutionSpaceOffsets);

    typedef typename boost::mpl::fold<
        typename boost::mpl::transform<
            /* This as_vector is probably not needed for boost::fusion 1.58
             * or higher. */
            typename result_of::as_vector<typename std::remove_reference<
                  decltype(bilinearForm_.getTerms())>::type
                >::type
          , mpl::second<boost::mpl::_1>
          >::type
      , boost::mpl::set0<>
      , boost::mpl::insert<boost::mpl::_1,boost::mpl::_2>
      >::type LFIndices;

    auto lfIndices = LFIndices{};
    for_each(lfIndices,
            localToGlobalRHSCopyHelper<decltype(solutionZip),
                                        decltype(cp)>
                                        (solutionZip, cp));

  }
}


template<class InnerProduct, class BilinearForm>
template <size_t spaceIndex, class ValueType>
void DPGSystemAssembler<InnerProduct, BilinearForm>::
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


template<class InnerProduct, class BilinearForm>
template <size_t spaceIndex, class ValueType>
void DPGSystemAssembler<InnerProduct, BilinearForm>::
applyDirichletBoundaryToMatrix
                      (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& boundaryValue)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  const size_t spaceSize =
        std::get<spaceIndex>(solutionSpaces_).size();

  size_t globalOffset;
  {
    /* set up global offsets */
    size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

    fold(zip(globalSolutionSpaceOffsets, solutionSpaces_),
         (size_t)0, globalOffsetHelper());

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


template<class InnerProduct, class BilinearForm>
template <size_t spaceIndex, class ValueType>
void DPGSystemAssembler<InnerProduct, BilinearForm>::
applyDirichletBoundaryToRhs
                      (BlockVector<FieldVector<double,1> >& rhs,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& boundaryValue)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  const size_t spaceSize =
        std::get<spaceIndex>(solutionSpaces_).size();

  size_t globalOffset;
  {
    /* set up global offsets */
    size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

    fold(zip(globalSolutionSpaceOffsets, solutionSpaces_),
         (size_t)0, globalOffsetHelper());

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

template<class InnerProduct, class BilinearForm>
template <size_t spaceIndex, unsigned int dim>
void DPGSystemAssembler<InnerProduct, BilinearForm>::
applyWeakBoundaryCondition
                    (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                     FieldVector<double, dim> beta,
                     double mu)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  typedef typename std::tuple_element<0,SolutionSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(solutionSpaces_).gridView();

  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  fold(zip(globalSolutionSpaceOffsets, solutionSpaces_),
       (size_t)0, globalOffsetHelper());
  size_t globalOffset = globalSolutionSpaceOffsets[spaceIndex];

  auto localView     = std::get<spaceIndex>(solutionSpaces_).localView();
  auto localIndexSet = std::get<spaceIndex>(solutionSpaces_).localIndexSet();

  for(const auto& e : elements(gridView))
  {
    localView.bind(e);
    localIndexSet.bind(localView);

    const auto& localFiniteElement = localView.tree().finiteElement();

    const unsigned int quadratureOrder
        = 2*localFiniteElement.localBasis().order();

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

            const double integrationWeight
              = intersection.geometry().integrationElement(quadFacePos)
              * quadFace[pt].weight();

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
                           * integrationWeight;
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


template<class InnerProduct, class BilinearForm>
template<size_t spaceIndex, unsigned int dim>
void DPGSystemAssembler<InnerProduct, BilinearForm>::
defineCharacteristicFaces(BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                          BlockVector<FieldVector<double,1> >& rhs,
                          const FieldVector<double,dim>& beta,
                          double delta)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  typedef typename std::tuple_element<spaceIndex,SolutionSpaces>::type::GridView GridView;
  GridView gridView = std::get<spaceIndex>(solutionSpaces_).gridView();

  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  fold(zip(globalSolutionSpaceOffsets, solutionSpaces_),
       (size_t)0, globalOffsetHelper());
  size_t globalOffset = globalSolutionSpaceOffsets[spaceIndex];

  auto solutionLocalView = at_c<spaceIndex>(solutionSpaces_).localView();
  auto localIndexSet = at_c<spaceIndex>(solutionSpaces_).localIndexSet();

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


template<class InnerProduct, class BilinearForm>
template <size_t spaceIndex, class MinInnerProduct, unsigned int dim>
void DPGSystemAssembler<InnerProduct, BilinearForm>::
applyMinimization
            (BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
             MinInnerProduct minInnerProduct,
             FieldVector<double, dim> beta,
             double delta,
             double epsilon)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  typedef typename std::tuple_element<spaceIndex,SolutionSpaces>::type::GridView GridView;
  GridView gridView = std::get<spaceIndex>(solutionSpaces_).gridView();

  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  fold(zip(globalSolutionSpaceOffsets, solutionSpaces_),
       (size_t)0, globalOffsetHelper());
  size_t globalOffset = globalSolutionSpaceOffsets[spaceIndex];

  size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

  // get local view for solution space (necessary if we want to use inner product) /TODO inefficient (why?)
  auto solutionLocalView = as_vector(transform(solutionSpaces_, getLocalView()));

  auto localIndexSet = at_c<spaceIndex>(solutionSpaces_).localIndexSet();

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

#endif // DUNE_DPG_SYSTEM_ASSEMBLER_DPG_HH
