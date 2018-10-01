// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_ERROR_TOOLS
#define DUNE_DPG_ERROR_TOOLS

#include <algorithm>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <vector>

#include <dune/common/hybridutilities.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dune/dpg/functions/localindexsetiteration.hh>
#include "assemble_helper.hh"
#include "quadrature.hh"

#if DUNE_DPG_USE_LEAST_SQUARES_INSTEAD_OF_CHOLESKY
#  include <dune/dpg/leastsquares.hh>
#else
#  include <dune/dpg/cholesky.hh>
#endif


namespace Dune {

  namespace detail {
    template<class GlobalVectorType, class LocalVectorType,
            class LocalIndexSets,
            class Offsets>
    inline void getLocalCoefficients(
        const GlobalVectorType& solution,
        LocalVectorType& solutionElement,
        const LocalIndexSets& localIndexSets,
        const Offsets& localOffsets,
        const Offsets& globalOffsets) {
      Hybrid::forEach(
          Std::make_index_sequence<
              std::tuple_size<LocalIndexSets>::value>{},
          [&](auto i) {
            auto const & localIndexSet = std::get<i>(localIndexSets);
            iterateOverLocalIndices(
              localIndexSet,
              [&](size_t j, auto gj) {
                solutionElement[j + localOffsets[i]]
                    = solution[globalOffsets[i] + gj[0]];
              },
              [&](size_t j){
                solutionElement[j + localOffsets[i]] = 0;
              },
              [&](size_t j, auto gj, double wj) {
                solutionElement[j + localOffsets[i]]
                    += wj * solution[globalOffsets[i] + gj[0]];
              }
            );
          });
    }
  }


  //*******************************************************************
  class ErrorTools
  {
  public:
    ErrorTools() = delete;

    template <class FEBasis>
    static double l2norm(
        const FEBasis& ,
        const BlockVector<FieldVector<double,1> >& );

    template <unsigned int subsamples, class FEBasis, class Function>
    static
    double computeL2error(const FEBasis& ,
                          const BlockVector<FieldVector<double,1> >& ,
                          Function&&,
                          unsigned int = 5);

    template <class InnerProduct, class LinearForm,
              class RhsFunction, class VectorType>
    static
    double aPosterioriL2Error(InnerProduct& ,
                              LinearForm& ,
                              const RhsFunction& ,
                              const VectorType&);

    template <class BilinearForm, class InnerProduct, class VectorType>
    static
    double aPosterioriError(BilinearForm& ,
                            InnerProduct& ,
                            const VectorType& ,
                            const VectorType& );

    template <class Grid, class EntitySeed>
    static
    double DoerflerMarking(Grid& ,
                           double ,
                           std::vector<std::tuple<EntitySeed, double>>&& );

    template <class BilinearForm, class InnerProduct, class VectorType>
    static
    std::vector<std::tuple<typename
        std::tuple_element_t<0,typename BilinearForm::SolutionSpaces>
          ::GridView::template Codim<0>::Entity::EntitySeed,
        double>>
    squaredCellwiseResidual(
        BilinearForm& ,
        InnerProduct& ,
        const VectorType& ,
        const VectorType& );

    template <class BilinearForm, class InnerProduct,
              class APosterioriInnerProduct, class LinearForm,
              class RhsFunction, class VectorType>
    static
    std::vector<std::tuple<typename
        std::tuple_element_t<0,typename BilinearForm::SolutionSpaces>
          ::GridView::template Codim<0>::Entity::EntitySeed,
        double>>
    squaredCellwiseResidual(
             BilinearForm& ,
             InnerProduct& ,
             APosterioriInnerProduct& ,
             LinearForm& ,
             const RhsFunction& ,
             const VectorType& ,
             const VectorType& ,
             double );

  private:
    template <class LocalView>
    static double l2normSquaredElement(
        const LocalView& ,
        const BlockVector<FieldVector<double,1> >& );

    template <unsigned int subsamples, class LocalView, class Function>
    static double computeL2errorSquareElement(
        const LocalView& ,
        const BlockVector<FieldVector<double,1> >& ,
        Function&&,
        unsigned int = 5);

    template <class InnerProduct, class LinearForm,
              class RhsFunction, class SolutionLocalViews,
              class VectorType>
    static
    double aPosterioriL2ErrorSquareElement(InnerProduct& ,
                                           LinearForm& ,
                                           const RhsFunction& ,
                                           SolutionLocalViews& ,
                                           VectorType&);

    template <class BilinearForm,class InnerProduct,
              class TestLocalViews, class SolutionLocalViews,
              class VectorType>
    static
    double aPosterioriErrorSquareElement(BilinearForm& ,
                                         InnerProduct& ,
                                         TestLocalViews& ,
                                         SolutionLocalViews& ,
                                         VectorType& ,
                                         VectorType& );
  };

//*******************************************************************

/**
 * \brief computes the square of the L2 of a given local FEfunction u
 *        over a given element.
 *
 * \param localView     the local view of the element
 * \param u             the vector containing the local FE coefficients
 */
  template <class LocalView>
  double ErrorTools::l2normSquaredElement(
      const LocalView& localView,
      const BlockVector<FieldVector<double,1> >& u)
  {
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    constexpr int dim = Element::mydimension;
    const auto geometry = element.geometry();

    auto&& localBasis = localView.tree().finiteElement().localBasis();

    // take 2*order since we integrate u^2
    const unsigned int quadratureOrder
        = 2*localBasis.order();

    const QuadratureRule<double, dim>& quad =
        QuadratureRules<double, dim>::rule(element.type(), quadratureOrder);

    double l2NormSquared = 0;

    for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();

      // The multiplicative factor in the integral transformation formula
      const double integrationElement = geometry.integrationElement(quadPos);

      std::vector<FieldVector<double,1> > shapeFunctionValues;
      localBasis.evaluateFunction(quadPos, shapeFunctionValues);

      double uQuad = 0;
      for(unsigned int i=0, imax=shapeFunctionValues.size(); i<imax; i++)
      {
        uQuad += shapeFunctionValues[i]*u[i];
      }

      l2NormSquared += uQuad * uQuad * quad[pt].weight() * integrationElement;
    }

    return l2NormSquared;
  }

/**
 * \brief Compute the L2 error of a FE function
 *
 * \param feBasis       the finite element basis
 * \param u             the vector containing the FE coefficients
 */
  template <class FEBasis>
  double ErrorTools::l2norm(const FEBasis& feBasis,
                            const BlockVector<FieldVector<double,1> >& u)
  {
    static_assert(!is_RefinedFiniteElement<FEBasis>::value,
        "l2norm only implemented for unrefined FE spaces!");

    auto gridView = feBasis.gridView();

    double l2NormSquared = 0.;

    auto localView = feBasis.localView();

    for(const auto& e : elements(gridView))
    {
      localView.bind(e);

      const size_t dofFEelement = localView.size();

      // We take the coefficients of u that correspond to the current
      // element e and store them in uElement.
      BlockVector<FieldVector<double,1> > uElement(dofFEelement);
      for (size_t i=0; i<dofFEelement; i++)
      {
          uElement[i] = u[ localView.index(i)[0] ];
      }

      l2NormSquared += l2normSquaredElement(localView, uElement);
    }

    return std::sqrt(l2NormSquared);
  }

/**
 * \brief Returns the computation in a given element of the L2 error
          between the exact solution uRef and the fem solution u.
 *
 * \tparam subsamples   number of subsamples (per direction) used for
 *                      the quadrature (must be >= 1)
 * \param localView     the local view of the element
 * \param u             the vector containing the computed solution
 * \param uRef          the expression for the exact solution.
 * \param quadOrder     polynomial degree up to which (x2) the quadrature
 *                      shall be exact.
 *                      If quadOrder < polynomial degree of the local finite
 *                      element, the polynomial degree of the local finite
 *                      element (x2) will be used instead.
 */
  template <unsigned int subsamples, class LocalView, class Function>
  double ErrorTools::computeL2errorSquareElement(
      const LocalView& localView,
      const BlockVector<FieldVector<double,1> >& u,
      Function&& uRef,
      const unsigned int quadOrder)
  {

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    constexpr int dim = Element::mydimension;
    const auto geometry = element.geometry();

    // Get set of shape functions for this element
    auto&& localBasis = localView.tree().finiteElement().localBasis();

    // taking twice the quadrature order as we integrate over a square
    const unsigned int quadratureOrder
        = std::max(2*quadOrder, 2*localBasis.order());

    const QuadratureRule<double, dim>& quadSection =
        QuadratureRules<double, dim>::rule(element.type(), quadratureOrder);
    const SubsampledQuadratureRule<double, subsamples, dim>& quad(quadSection);

    // Variables employed in the loop
    double errSquare = 0;     // we store here the square of the error
    double uQuad = 0;         // we store here the value of u at a quadrature point

    // Loop over all quadrature points
    for (size_t pt=0; pt < quad.size(); pt++) {

      uQuad = 0;

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();
      // Position of the current quadrature point in the current element
      const FieldVector<double,dim> globalQuadPos = geometry.global(quadPos);

      // The multiplicative factor in the integral transformation formula
      const double integrationElement = geometry.integrationElement(quadPos);

      // Evaluate all shape function values at quadPos (which is a
      // quadrature point in the reference element)
      std::vector<FieldVector<double,1> > shapeFunctionValues;
      localBasis.evaluateFunction(quadPos, shapeFunctionValues);

      // Evaluation of u at the point globalQuadPos, which is quadPos
      // mapped to the physical domain
      for(unsigned int i=0; i<shapeFunctionValues.size(); i++)
      {
        uQuad += shapeFunctionValues[i]*u[i];
      }

      // Value of uRef at globalQuadPos
      const double uExactQuad = uRef(globalQuadPos);

      // we add the squared error at the quadrature point
      errSquare += (uQuad - uExactQuad)*(uQuad - uExactQuad)
                 * quad[pt].weight() * integrationElement;
    }

    return errSquare;
  }

/**
 * \brief Computation in the whole mesh of the L2 error
          between the exact solution uRef and the fem solution u.
 *
 * \tparam subsamples   number of subsamples (per direction) used for
 *                      the quadrature
 *                      (must be >= 1)
 * \param feBasis       the finite element basis
 * \param u             the vector containing the computed solution
 * \param uRef          the expression for the exact solution.
 * \param quadratureOrder
 *                      polynomial degree up to which (x2) the quadrature
 *                      shall be exact.
 *                      If quadOrder < polynomial degree of the local finite
 *                      element, the polynomial degree of the local finite
 *                      element (x2) will be used instead.
 */
  template <unsigned int subsamples, class FEBasis, class Function>
  double ErrorTools::computeL2error(
      const FEBasis& feBasis,
      const BlockVector<FieldVector<double,1> >& u,
      Function&& uRef,
      const unsigned int quadratureOrder)
  {
    // Get the grid view from the finite element basis
    typedef typename FEBasis::GridView GridView;
    GridView gridView = feBasis.gridView();

    // Variables where we will store the errors
    double errSquare = 0.;

    // A view on the FE basis on a single element
    auto localView = feBasis.localView();

    // A loop over all elements of the grid
    for(const auto& e : elements(gridView))
    {
      // Bind the local FE basis view to the current element
      localView.bind(e);

      // Now we take the coefficients of u that correspond to the
      // current element e. They are stored in uElement.
      // number of dofs of the finite element inside the element
      // (remark: this value will vary if we do p-refinement)
      size_t dofFEelement = localView.size();

      // We take the coefficients of u that correspond to the current
      // element e. They are stored in uElement.
      BlockVector<FieldVector<double,1> > uElement(dofFEelement);
      for (size_t i=0; i<dofFEelement; i++)
      {
          uElement[i] = u[ localView.index(i)[0] ];
      }
      // Now we compute the error inside the element
      errSquare += computeL2errorSquareElement<subsamples>
                       (localView, uElement, uRef, quadratureOrder);
    }

    return std::sqrt(errSquare);
  }

/**
 * \brief Computation of an alternative a posteriori error in (u,theta)
 *
 * for a detailed description see ErrorTools::aPosterioriL2Error
 *
 * \param innerProduct            the inner product
 * \param linearForm              the linear form
 * \param f                       the right-hand side function
 * \param solutionLocalViews      local views for the solution
 * \param solution                the computed solution
 */
  template <class InnerProduct, class LinearForm,
            class RhsFunction, class SolutionLocalViews,
            class VectorType>
  double ErrorTools::aPosterioriL2ErrorSquareElement(
      InnerProduct& innerProduct,
      LinearForm& linearForm,
      const RhsFunction& f,
      SolutionLocalViews& solutionLocalViews,
      VectorType& solution)
  {
    using namespace Dune::detail;

    typedef typename InnerProduct::TestSpaces SolutionSpaces;

    // Create and fill vector with offsets for global dofs
    size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
    computeOffsets(globalSolutionSpaceOffsets, *innerProduct.getTestSpaces());

    // Create and fill vector with offsets for local dofs on element
    size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
    const size_t localSolutionDofs
        = computeOffsets(localSolutionSpaceOffsets, solutionLocalViews);

    // Create and fill vectors with cofficients
    // corresponding to local dofs on element
    // for the solution and for the righthand side
    BlockVector<FieldVector<double,1> > solutionElement(localSolutionDofs);

    detail::getLocalCoefficients(solution, solutionElement,
        solutionLocalViews,
        localSolutionSpaceOffsets, globalSolutionSpaceOffsets);

    double errSquare = 0;

    // We grab the inner product matrix in the innerProductMatrix variable (IP)
    Matrix<FieldMatrix<double,1,1> > innerProductMatrix;
    innerProduct.bind(solutionLocalViews);
    innerProduct.getLocalMatrix(innerProductMatrix);

    // We compute solutionElement^T * IP * solutionElement
    BlockVector<FieldVector<double,1> > tmpVector(solutionElement);
    tmpVector = 0;
    innerProductMatrix.mv(solutionElement, tmpVector);
    errSquare += tmpVector * solutionElement;

    // We grab the linear term depending on B and on the rhs f
    BlockVector<FieldVector<double,1> > linearFormVector;
    linearForm.bind(solutionLocalViews);
    linearForm.getLocalVector(linearFormVector);
    errSquare += linearFormVector * solutionElement;

    // We grab the term depending only on the rhs f
    // TODO Not tested for non-constant f !!!!!!!!
    double tmpValue = 0;
    // TODO adjust quadrature for non-constant RHS
    unsigned int quadratureOrder = 5;
    auto element = std::get<0>(solutionLocalViews).element();
    const auto geometry = element.geometry();
    typedef decltype(element) Element;
    constexpr int dim = Element::mydimension;
    const Dune::QuadratureRule<double, dim>& quad =
            Dune::QuadratureRules<double, dim>::rule(element.type(),
                                                     quadratureOrder);
    for ( size_t pt=0, qsize=quad.size(); pt < qsize; pt++ ) {
      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();
      // Global position of the current quadrature point
      const FieldVector<double,dim> globalQuadPos
          = geometry.global(quadPos);
      tmpValue+= geometry.integrationElement(quadPos)
                 * quad[pt].weight()
                 * f(globalQuadPos)
                 * f(globalQuadPos);
    }
    errSquare += tmpValue;

    return errSquare;
  }

/**
 * \brief Computation of an alternative a posteriori error in (u,theta)
 *
 * Computation of an alternative a posteriori error of the form
 * (< solution, solution >_{innerProduct} + linearForm(solution) + f^2)^{1/2}
 * in the transport case, this is used to compute
 * ||u-w||_L_2^2 + ||Bw-f||_L_2^2
 * where w is the lifting of theta (so we have to use a formulation
 * where we use the lifting and not theta itself)
 * and B is the operator of the conforming formulation
 *
 * In this example, we have
 * innerProduct = < u-w,u-w>_L_2 + < Bw, BW>_L_2
 * linearForm = -2 < Bw, f>_L_2
 * constant term f = f
 *
 * \param innerProduct       the inner product
 * \param linearForm         the linear form
 * \param f                  the right-hand side function
 * \param solution           the computed solution
 */
  template <class InnerProduct, class LinearForm,
            class RhsFunction, class VectorType>
  double ErrorTools::aPosterioriL2Error(InnerProduct& innerProduct,
                                        LinearForm& linearForm,
                                        const RhsFunction& f,
                                        const VectorType& solution)
  {
    using namespace Dune::detail;

    typedef typename std::tuple_element<0, typename InnerProduct::TestSpaces>
                        ::type::GridView GridView;
    const GridView gridView
        = std::get<0>(*innerProduct.getTestSpaces()).gridView();

    typedef typename InnerProduct::TestSpaces SolutionSpaces;

    typedef getLocalViews_t<SolutionSpaces>  SolutionLocalViews;

    SolutionLocalViews solutionLocalViews
        = getLocalViews(*innerProduct.getTestSpaces());

    double squaredResidual = 0.;

    for(const auto& e : elements(gridView))
    {
      bindLocalViews(solutionLocalViews, e);

      squaredResidual
          += aPosterioriL2ErrorSquareElement(innerProduct,
                                             linearForm,
                                             f,
                                             solutionLocalViews,
                                             solution);
    }

    return std::sqrt(squaredResidual);
  }



/**
 * \brief Computation of a posteriori error in (u, theta)
 *
 * \param bilinearForm        the bilinear form
 * \param innerProduct        the inner product
 * \param testLocalViews
 * \param solutionLocalViews
 * \param testLocalIndexSets
 * \param solutionLocalIndexSets
 * \param solution            the computed solution
 * \param rhs                 the right-hand side
 */
  template <class BilinearForm,class InnerProduct,
            class TestLocalViews, class SolutionLocalViews,
            class VectorType>
  double ErrorTools::aPosterioriErrorSquareElement(
      BilinearForm& bilinearForm,
      InnerProduct& innerProduct,
      TestLocalViews& testLocalViews,
      SolutionLocalViews& solutionLocalViews,
      VectorType& solution,
      VectorType& rhs)
  {
    using namespace Dune::detail;

    typedef typename BilinearForm::SolutionSpaces SolutionSpaces;
    typedef typename BilinearForm::TestSpaces EnrichedTestspaces;

    // Create and fill vector with offsets for global dofs
    size_t globalTestSpaceOffsets[std::tuple_size<EnrichedTestspaces>::value];
    size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

    computeOffsets(globalTestSpaceOffsets, *bilinearForm.getTestSpaces());
    computeOffsets(globalSolutionSpaceOffsets,
                   *bilinearForm.getSolutionSpaces());

    // Create and fill vector with offsets for local dofs on element
    size_t localTestSpaceOffsets[std::tuple_size<EnrichedTestspaces>::value];
    size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

    const size_t localSolutionDofs
        = computeOffsets(localSolutionSpaceOffsets, solutionLocalViews);
    const size_t localTestDofs
        = computeOffsets(localTestSpaceOffsets, testLocalViews);

    // Create and fill vectors with cofficients
    // corresponding to local dofs on element
    // for the solution and for the righthand side
    BlockVector<FieldVector<double,1> > solutionElement(localSolutionDofs);
    BlockVector<FieldVector<double,1> > rhsElement(localTestDofs);

    detail::getLocalCoefficients(solution, solutionElement,
        solutionLocalViews,
        localSolutionSpaceOffsets, globalSolutionSpaceOffsets);
    detail::getLocalCoefficients(rhs, rhsElement,
        testLocalViews,
        localTestSpaceOffsets, globalTestSpaceOffsets);

    // We grab the inner product matrix in the innerProductMatrix variable (IP)
    Matrix<FieldMatrix<double,1,1> > innerProductMatrix;
    innerProduct.bind(testLocalViews);
    innerProduct.getLocalMatrix(innerProductMatrix);

    // We grab the bilinear form matrix in the bilinearFormMatrix variable
    Matrix<FieldMatrix<double,1,1> > bilinearFormMatrix;
    bilinearForm.bind(testLocalViews, solutionLocalViews);
    bilinearForm.getLocalMatrix(bilinearFormMatrix);

    // compute Bu - f (we do f-= Bu so the output is in rhsElement)
    bilinearFormMatrix.mmv(solutionElement,rhsElement);

    // Solve for the Riesz lift and then compute the residual
    BlockVector<FieldVector<double,1> > tmpVector(rhsElement);
#if DUNE_DPG_USE_LEAST_SQUARES_INSTEAD_OF_CHOLESKY
    solveLeastSquares(innerProductMatrix, tmpVector);
#else
    {
      Cholesky<Matrix<FieldMatrix<double,1,1>>> cholesky(innerProductMatrix);
      cholesky.apply(tmpVector);
    }
#endif

    // compute the scalar product of the following two vectors
    return tmpVector * rhsElement;
  }

/**
 * \brief Computation of a posteriori error in (u, theta)
 *
 * \param bilinearForm        the bilinear form
 * \param innerProduct        the inner product
 * \param solution            the computed solution
 * \param rhs                 the right-hand side
 */
  template <class BilinearForm,class InnerProduct,
            class VectorType>
  double ErrorTools::aPosterioriError(BilinearForm& bilinearForm,
                                      InnerProduct& innerProduct,
                                      const VectorType& solution,
                                      const VectorType& rhs)
  {
    using namespace Dune::detail;

    typedef typename std::tuple_element
        < 0
        , typename BilinearForm::SolutionSpaces
        >::type::GridView GridView;
    const GridView gridView
        = std::get<0>(*bilinearForm.getSolutionSpaces()).gridView();

    typedef typename BilinearForm::SolutionSpaces SolutionSpaces;
    typedef typename BilinearForm::TestSpaces EnrichedTestspaces;

    typedef detail::getLocalViews_t<SolutionSpaces>  SolutionLocalViews;
    typedef detail::getLocalViews_t<EnrichedTestspaces>  TestLocalViews;

    SolutionLocalViews solutionLocalViews
        = getLocalViews(*bilinearForm.getSolutionSpaces());
    TestLocalViews testLocalViews
        = getLocalViews(*bilinearForm.getTestSpaces());

    double squaredResidual = 0.;

    for(const auto& e : elements(gridView))
    {
      bindLocalViews(testLocalViews, e);
      bindLocalViews(solutionLocalViews, e);

      squaredResidual
          += aPosterioriErrorSquareElement(bilinearForm,
                                           innerProduct,
                                           testLocalViews,
                                           solutionLocalViews,
                                           solution,
                                           rhs);
   }

   return std::sqrt(squaredResidual);

  }

/**
 * \brief Mark elements for refinement according to Dörfler's strategy
 *
 * This means, given a ratio \f$\theta \in (0,1]\f$, the set
 * \f$\mathcal M \subset \Omega_h\f$ of marked elements satisfies
 * \f[\mathrm{err}(u_h, \mathcal M)
     \geq \mathrm{err}(u_h, \Omega_h).\f]
 *
 * \param grid
 * \param ratio  the marking ratio \f$\theta \in (0,1]\f$
 * \param errorEstimates holds a squared error estimate for each entity
 *
 * \return estimate for the squared global a posteriori error
 */
  template <class Grid, class EntitySeed>
  double ErrorTools::DoerflerMarking(
      Grid& grid,
      double ratio,
      std::vector<std::tuple<EntitySeed, double>>&& errorEstimates)
  {
    using namespace Dune::detail;

    using Entity = typename Grid::template Codim<0>::Entity;
    static_assert(std::is_same<EntitySeed,
        typename Entity::EntitySeed>::value,
        "EntitySeed type does not fit the Grid type.");

    std::sort(errorEstimates.begin(), errorEstimates.end(),
        [](const auto& a, const auto& b) {
            return std::get<1>(a) > std::get<1>(b);
        });
    const double errorSquared = std::accumulate(
        errorEstimates.cbegin(), errorEstimates.cend(), 0.,
        [](double a, const auto& b) {
            return a + std::get<1>(b);
        });
    // Since we used squared errors, we have to square the marking ratio
    const double targetError = ratio * ratio * errorSquared;
    double error = 0.;
    auto currElem = errorEstimates.begin();
    while(error < targetError) {
      error += std::get<1>(*currElem);
      grid.mark(1, grid.entity(std::get<0>(*currElem)));
      ++currElem;
    }

    return errorSquared;
  }

/**
 * \brief Element-wise squared residual of an (ultra-)weak formulation
 *
 * \todo document \p aPosterioriInnerProduct, \p linearForm and \p splitRatio
 *
 * \param bilinearForm
 * \param innerProduct
 * \param aPosterioriInnerProduct
 * \param linearForm
 * \param f
 * \param solution  FE solution
 * \param rhs  Rhs of the problem in enriched test space
 * \param splitRatio
 *
 * \return cell-wise estimate for the squared a posteriori error
 */
  template <class BilinearForm, class InnerProduct,
            class APosterioriInnerProduct, class LinearForm,
            class RhsFunction, class VectorType>
  std::vector<std::tuple<typename
      std::tuple_element_t<0,typename BilinearForm::SolutionSpaces>
        ::GridView::template Codim<0>::Entity::EntitySeed,
      double>>
  ErrorTools::squaredCellwiseResidual(
      BilinearForm& bilinearForm,
      InnerProduct& innerProduct,
      APosterioriInnerProduct& aPosterioriInnerProduct,
      LinearForm& linearForm,
      const RhsFunction& f,
      const VectorType& solution,
      const VectorType& rhs,
      double splitRatio)
  {
    using namespace Dune::detail;

    assert(splitRatio >= 0 && splitRatio <= 1);
    using GridView
        = typename std::tuple_element<0,typename BilinearForm::SolutionSpaces>
                        ::type::GridView;
    using Entity = typename GridView::template Codim<0>::Entity;
    using EntitySeed = typename Entity::EntitySeed;

    const GridView gridView
        = std::get<0>(*bilinearForm.getSolutionSpaces()).gridView();

    using SolutionSpaces = typename BilinearForm::SolutionSpaces;
    using EnrichedTestspaces = typename BilinearForm::TestSpaces;

    using SolutionLocalViews = getLocalViews_t<SolutionSpaces>;
    using TestLocalViews = getLocalViews_t<EnrichedTestspaces>;

    SolutionLocalViews solutionLocalViews
        = getLocalViews(*bilinearForm.getSolutionSpaces());
    TestLocalViews testLocalViews
        = getLocalViews(*bilinearForm.getTestSpaces());

    std::vector<std::tuple<EntitySeed, double>> errorEstimates;
    errorEstimates.reserve(gridView.size(0));
    for(const auto& e : elements(gridView))
    {
      bindLocalViews(testLocalViews, e);
      bindLocalViews(solutionLocalViews, e);

      // Now we compute the error inside the element
      double elementError = 0;
      if (splitRatio > 0)
      {
        elementError += splitRatio*aPosterioriErrorSquareElement(
                                bilinearForm,
                                innerProduct,
                                testLocalViews,
                                solutionLocalViews,
                                solution,
                                rhs);
      }
      if (splitRatio < 1)
      {
        elementError += (1-splitRatio)*aPosterioriL2ErrorSquareElement(
                                aPosterioriInnerProduct,
                                linearForm,
                                f,
                                solutionLocalViews,
                                solution);
      }
      errorEstimates.emplace_back(e.seed(), elementError);
    }
    return errorEstimates;
  }

/**
 * \brief Element-wise squared residual of an (ultra-)weak formulation
 *
 * \param bilinearForm
 * \param innerProduct
 * \param solution  FE solution
 * \param rhs  Rhs of the problem in enriched test space
 *
 * \return cell-wise estimate for the squared a posteriori error
 */
  template <class BilinearForm, class InnerProduct, class VectorType>
  std::vector<std::tuple<typename
      std::tuple_element_t<0,typename BilinearForm::SolutionSpaces>
        ::GridView::template Codim<0>::Entity::EntitySeed,
      double>>
  ErrorTools::squaredCellwiseResidual(
      BilinearForm& bilinearForm,
      InnerProduct& innerProduct,
      const VectorType& solution,
      const VectorType& rhs)
  {
    using namespace Dune::detail;

    using GridView
        = typename std::tuple_element<0,typename BilinearForm::SolutionSpaces>
                        ::type::GridView;
    using Entity = typename GridView::template Codim<0>::Entity;
    using EntitySeed = typename Entity::EntitySeed;

    const GridView gridView
        = std::get<0>(*bilinearForm.getSolutionSpaces()).gridView();

    using SolutionSpaces = typename BilinearForm::SolutionSpaces;
    using EnrichedTestspaces = typename BilinearForm::TestSpaces;

    using SolutionLocalViews = getLocalViews_t<SolutionSpaces>;
    using TestLocalViews = getLocalViews_t<EnrichedTestspaces>;

    SolutionLocalViews solutionLocalViews
        = getLocalViews(*bilinearForm.getSolutionSpaces());
    TestLocalViews testLocalViews
        = getLocalViews(*bilinearForm.getTestSpaces());

    std::vector<std::tuple<EntitySeed, double>> errorEstimates;
    errorEstimates.reserve(gridView.size(0));
    for(const auto& e : elements(gridView))
    {
      bindLocalViews(testLocalViews, e);
      bindLocalViews(solutionLocalViews, e);

      // Now we compute the error inside the element
      double elementError = aPosterioriErrorSquareElement(
                                bilinearForm,
                                innerProduct,
                                testLocalViews,
                                solutionLocalViews,
                                solution,
                                rhs);
      errorEstimates.emplace_back(e.seed(), elementError);
    }
    return errorEstimates;
  }

} // end namespace Dune

#endif // DUNE_DPG_ERROR_TOOLS
