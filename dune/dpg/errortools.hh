// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_ERROR_TOOLS
#define DUNE_DPG_ERROR_TOOLS

#include <algorithm>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <vector>

#include <boost/fusion/algorithm/transformation/transform.hpp>
#include <boost/fusion/algorithm/transformation/zip.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/functional/generation/make_fused_procedure.hpp>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include "assemble_helper.hh"
#include "cholesky.hh"
#include "quadrature.hh"


namespace Dune {


  //*******************************************************************
  class ErrorTools
  {
  public:
    ErrorTools() {};
    template <unsigned int subsamples, class LocalView, class VolumeTerms>
    double computeL2errorSquareElement(const LocalView& ,
                                 BlockVector<FieldVector<double,1> >& ,
                                 VolumeTerms&&,
                                 unsigned int = 5);

    template <unsigned int subsamples, class FEBasis, class VolumeTerms>
    double computeL2error(const FEBasis& ,
                          BlockVector<FieldVector<double,1> >& ,
                          VolumeTerms&&,
                          unsigned int = 5);

    template <class InnerProduct, class LinearForm,
              class RhsFunction, class LocalViewsSolution,
              class LocalIndexSetSolution, class VectorType>
    double aPosterioriL2ErrorSquareElement(InnerProduct& ,
                                      LinearForm& ,
                                      const RhsFunction& ,
                                      LocalViewsSolution& ,
                                      LocalIndexSetSolution& ,
                                      VectorType&);

    template <class InnerProduct, class LinearForm,
              class RhsFunction, class VectorType>
    double aPosterioriL2Error(InnerProduct& ,
                                        LinearForm& ,
                                        const RhsFunction& ,
                                        const VectorType&);

    template <class BilinearForm,class InnerProduct,
              class LocalViewsTest, class LocalViewsSolution,
              class LocalIndexSetTest, class LocalIndexSetSolution,
              class VectorType>
    double aPosterioriErrorSquareElement(BilinearForm& ,
                                      InnerProduct& ,
                                      LocalViewsTest& ,
                                      LocalViewsSolution& ,
                                      LocalIndexSetTest& ,
                                      LocalIndexSetSolution& ,
                                      VectorType& ,
                                      VectorType& );

    template <class BilinearForm, class InnerProduct, class VectorType>
    double aPosterioriError(BilinearForm& ,
                            InnerProduct& ,
                            const VectorType& ,
                            const VectorType& );

    template <class Grid, class BilinearForm, class InnerProduct,
              class APosterioriInnerProduct, class LinearForm,
              class RhsFunction, class VectorType>
    double DoerflerMarking(Grid& ,
                           double ,
                           BilinearForm& ,
                           InnerProduct& ,
                           APosterioriInnerProduct& ,
                           LinearForm& ,
                           const RhsFunction& ,
                           const VectorType& ,
                           const VectorType& ,
                           double );
  };

//*******************************************************************

/**
 * \brief Returns the computation in a given element of the L2 error
          between the exact solution uRef and the fem solution u.
 *
 * \tparam subsamples   number of subsamples (per direction) used for the quadrature
 *                      (must be >= 1)
 * \param localView      the local view of the element
 * \param u              the vector containing the computed solution
 * \param uRef           the expression for the exact solution.
 * \param quadOrder      polynomial degree up to which (x2) the quadrature shall be exact.
 *                       If quadOrder < polynomial degree of the local finite element, the
 *                       polynomial degree of the local finite element (x2)
 *                       will be used instead
 */
  template <unsigned int subsamples, class LocalView,class VolumeTerms>
  double ErrorTools::computeL2errorSquareElement(const LocalView& localView,
                                           BlockVector<FieldVector<double,1> >& u,
                                           VolumeTerms&& uRef,
                                           unsigned int quadOrder
                                          )
  {

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    const int dim = Element::dimension;
    auto geometry = element.geometry();

    // Get set of shape functions for this element
    const auto& localFiniteElement = localView.tree().finiteElement();

    const unsigned int quadratureOrder = std::max(2*quadOrder, 2*localFiniteElement.localBasis().order()); // Remark: the quadrature order has to be an even number

    const QuadratureRule<double, dim>& quadSection =
        QuadratureRules<double, dim>::rule(element.type(), quadratureOrder);
    const SubsampledQuadratureRule<double, subsamples, dim>& quad(quadSection);

    // Variables employed in the loop
    double errSquare = 0;     // we store here the square of the error
    double uQuad = 0;         // we store here the value of u at a quadrature point

    // Loop over all quadrature points
    for (size_t pt=0; pt < quad.size(); pt++) {

      //Restart value uQuad
      uQuad = 0;

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();
      // Position of the current quadrature point in the current element
      const FieldVector<double,dim>& mapQuadPos = geometry.global(quadPos);  // we get the global coordinate of quadPos

      // The multiplicative factor in the integral transformation formula
      const double integrationElement = element.geometry().integrationElement(quadPos);

      // Evaluate all shape function values at quadPos (which is a quadrature point in the reference element)
      std::vector<FieldVector<double,1> > shapeFunctionValues;
      localFiniteElement.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

      // Evaluation of u at the point mapQuadPos, which is quadPos mapped to the physical domain
      for(int i=0; i<shapeFunctionValues.size(); i++)
      {
        uQuad += shapeFunctionValues[i]*u[i];
      }

      // Value of uRef at mapQuadPos
      double uExactQuad = std::get<0>(uRef)(mapQuadPos);

      // we add the squared error at the quadrature point
      errSquare += (uQuad - uExactQuad)*(uQuad - uExactQuad) * quad[pt].weight() * integrationElement;
    }

    return errSquare;
  }

/**
 * \brief Computation in the whole mesh of the L2 error
          between the exact solution uRef and the fem solution u.
 *
 * \tparam subsamples   number of subsamples (per direction) used for the quadrature
 *                      (must be >= 1)
 * \param feBasis        the finite element basis
 * \param u              the vector containing the computed solution
 * \param uRef           the expression for the exact solution.
 * \param quadOrder      polynomial degree up to which (x2) the quadrature shall be exact.
 *                       If quadOrder < polynomial degree of the local finite element, the
 *                       polynomial degree of the local finite element (x2)
 *                       will be used instead
 */
  template <unsigned int subsamples, class FEBasis,class VolumeTerms>
  double ErrorTools::computeL2error(const FEBasis& feBasis,
                                    BlockVector<FieldVector<double,1> >& u,
                                    VolumeTerms&& uRef,
                                    unsigned int quadratureOrder
                                   )
  {
    // Get the grid view from the finite element basis
    typedef typename FEBasis::GridView GridView;
    GridView gridView = feBasis.gridView();

    // Variables where we will store the errors
    double errSquare = 0.;

    // A view on the FE basis on a single element
    auto localView = feBasis.localView();
    auto localIndexSet = feBasis.localIndexSet();

    // A loop over all elements of the grid
    for(const auto& e : elements(gridView))
    {
      // Bind the local FE basis view to the current element
      localView.bind(e);
      localIndexSet.bind(localView);

      // We get the global index of the current element
      int indexElement = gridView.indexSet().index(e);

      // Now we take the coefficients of u that correspond to the current element e. They are stored in uElement.
      // dof of the finite element inside the element (remark: this value will vary if we do p-refinement)
      size_t dofFEelement = localView.size();

      // We take the coefficients of u that correspond to the current element e. They are stored in uElement.
      BlockVector<FieldVector<double,1> > uElement(dofFEelement);
      for (size_t i=0; i<dofFEelement; i++)
      {
          uElement[i] = u[ localIndexSet.index(i)[0] ];
      }
      // Now we compute the error inside the element
      errSquare += computeL2errorSquareElement<subsamples>
                       (localView,uElement,uRef,quadratureOrder);
    }

    return std::sqrt(errSquare);
  }

/**
 * \brief Computation of an alternative a posteriori error in (u,theta)
 *
 * for a detailed description see ErrorTools::aPosterioriL2Error
 *
 * \param innerProduct          the inner product
 * \param linearForm            the linear form
 * \param f                     the right-hand side function
 * \param localViewsSolution    local views for the solution
 * \param solutionLocalIndexSets local index sets for the solution
 * \param solution              the computed solution
 */
  template <class InnerProduct, class LinearForm,
            class RhsFunction, class LocalViewsSolution,
            class LocalIndexSetSolution, class VectorType>
  double ErrorTools::aPosterioriL2ErrorSquareElement(InnerProduct& innerProduct,
                                      LinearForm& linearForm,
                                      const RhsFunction& f,
                                      LocalViewsSolution& localViewsSolution,
                                      LocalIndexSetSolution& solutionLocalIndexSets,
                                      VectorType& solution)
  {
    using namespace boost::fusion;
    using namespace Dune::detail;

    typedef typename InnerProduct::TestSpaces SolutionSpaces;

    // Create and fill vector with offsets for global dofs
    size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

    fold(zip(globalSolutionSpaceOffsets, innerProduct.getTestSpaces()),
           (size_t)0, globalOffsetHelper());

    // Create and fill vector with offsets for local dofs on element
    size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

    size_t localSolutionDofs = fold(zip(localSolutionSpaceOffsets,
                                        localViewsSolution),
                                    (size_t)0, globalOffsetHelper());
    // Create and fill vectors with cofficients
    // corresponding to local dofs on element
    // for the solution and for the righthand side
    BlockVector<FieldVector<double,1> > solutionElement(localSolutionDofs);

    for_each(zip(localViewsSolution, solutionLocalIndexSets,
                 localSolutionSpaceOffsets, globalSolutionSpaceOffsets),
         getLocalCoefficients<VectorType, BlockVector<FieldVector<double,1> > >(solution, solutionElement));

    double errSquare = 0;

    // We grab the inner product matrix in the innerProductMatrix variable (IP)
    Matrix<FieldMatrix<double,1,1> > innerProductMatrix;
    innerProduct.bind(localViewsSolution);
    innerProduct.getLocalMatrix(innerProductMatrix);

    // We compute solutionElement^T * IP * solutionElement
    BlockVector<FieldVector<double,1> > tmpVector(solutionElement);
    tmpVector = 0;
    innerProductMatrix.mv(solutionElement, tmpVector);
    errSquare += tmpVector * solutionElement;

    // We grab the linear term depending on B and on the rhs f
    BlockVector<FieldVector<double,1> > linearFormVector;
    linearForm.bind(localViewsSolution);
    linearForm.getLocalVector(linearFormVector);
    errSquare += linearFormVector * solutionElement;

    // We grab the term depending only on the rhs f
    // TODO Not tested for non-constant f !!!!!!!!
    double tmpValue = 0;
    unsigned int quadratureOrder = 5;    //TODO adjust quadrature for non-constant RHS
    auto element = at_c<0>(localViewsSolution).element();
    auto geometry = element.geometry();
    typedef decltype(element) Element;
    const int dim = Element::dimension;
    const Dune::QuadratureRule<double, dim>& quad =
            Dune::QuadratureRules<double, dim>::rule(element.type(),
                                                     quadratureOrder);
    for ( size_t pt=0, qsize=quad.size(); pt < qsize; pt++ ) {
      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();
      // Global position of the current quadrature point
      const FieldVector<double,dim>& globalQuadPos
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
    using namespace boost::fusion;
    using namespace Dune::detail;

    typedef typename std::tuple_element<0,typename InnerProduct::TestSpaces>
                        ::type::GridView GridView;
    const GridView gridView = std::get<0>(innerProduct.getTestSpaces()).gridView();

    typedef typename InnerProduct::TestSpaces SolutionSpaces;

    typedef typename result_of::as_vector<typename result_of::transform<SolutionSpaces,getLocalView>::type>::type SolutionLocalViews;

    SolutionLocalViews localViewsSolution = as_vector(transform(innerProduct.getTestSpaces(),getLocalView()));


    // We get the local index sets of the solution spaces
    auto solutionLocalIndexSets
            = as_vector(transform(innerProduct.getTestSpaces(),
                                  getLocalIndexSet()));

    // Variable where we compute the residual
    double res = 0.;

    for(const auto& e : elements(gridView))
    {
      // Bind localViews and localIndexSets
      for_each(localViewsSolution, applyBind<decltype(e)>(e));
      for_each(zip(solutionLocalIndexSets, localViewsSolution),
               make_fused_procedure(bindLocalIndexSet()));

      res += aPosterioriL2ErrorSquareElement(innerProduct,
                                             linearForm,
                                             f,
                                             localViewsSolution,
                                             solutionLocalIndexSets,
                                             solution);
   }

   return std::sqrt(res);

  }



/**
 * \brief Computation of a posteriori error in (u,theta)
 *
 * \param bilinearForm        the bilinear form
 * \param innerProduct        the inner product
 * \param solution           the computed solution
 * \param rhs                 the right-hand side
 */
  template <class BilinearForm,class InnerProduct,
            class LocalViewsTest, class LocalViewsSolution,
            class LocalIndexSetTest, class LocalIndexSetSolution,
            class VectorType>
  double ErrorTools::aPosterioriErrorSquareElement(BilinearForm& bilinearForm,
                                      InnerProduct& innerProduct,
                                      LocalViewsTest& localViewsTest,
                                      LocalViewsSolution& localViewsSolution,
                                      LocalIndexSetTest& testLocalIndexSets,
                                      LocalIndexSetSolution& solutionLocalIndexSets,
                                      VectorType& solution,
                                      VectorType& rhs)
  {
    using namespace boost::fusion;
    using namespace Dune::detail;

    typedef typename BilinearForm::SolutionSpaces SolutionSpaces;
    typedef typename BilinearForm::TestSpaces EnrichedTestspaces;

    // Create and fill vector with offsets for global dofs
    size_t globalTestSpaceOffsets[std::tuple_size<EnrichedTestspaces>::value];
    size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

    fold(zip(globalTestSpaceOffsets, bilinearForm.getTestSpaces()),
           (size_t)0, globalOffsetHelper());
    fold(zip(globalSolutionSpaceOffsets, bilinearForm.getSolutionSpaces()),
           (size_t)0, globalOffsetHelper());

    // Create and fill vector with offsets for local dofs on element
    size_t localTestSpaceOffsets[std::tuple_size<EnrichedTestspaces>::value];
    size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];

    size_t localSolutionDofs = fold(zip(localSolutionSpaceOffsets,
                                        localViewsSolution),
                                    (size_t)0, globalOffsetHelper());
    size_t localTestDofs = fold(zip(localTestSpaceOffsets,
                                    localViewsTest),
                                (size_t)0, globalOffsetHelper());
    // Create and fill vectors with cofficients
    // corresponding to local dofs on element
    // for the solution and for the righthand side
    BlockVector<FieldVector<double,1> > solutionElement(localSolutionDofs);
    BlockVector<FieldVector<double,1> > rhsElement(localTestDofs);

    for_each(zip(localViewsSolution, solutionLocalIndexSets,
                 localSolutionSpaceOffsets, globalSolutionSpaceOffsets),
         getLocalCoefficients<VectorType, BlockVector<FieldVector<double,1> > >(solution, solutionElement));
    for_each(zip(localViewsTest, testLocalIndexSets,
                 localTestSpaceOffsets, globalTestSpaceOffsets),
         getLocalCoefficients<VectorType, BlockVector<FieldVector<double,1> > >(rhs, rhsElement));

    // We grab the inner product matrix in the innerProductMatrix variable (IP)
    Matrix<FieldMatrix<double,1,1> > innerProductMatrix;
    innerProduct.bind(localViewsTest);
    innerProduct.getLocalMatrix(innerProductMatrix);

    // We grab the bilinear form matrix in the bilinearFormMatrix variable
    Matrix<FieldMatrix<double,1,1> > bilinearFormMatrix;
    bilinearForm.bind(localViewsTest, localViewsSolution);
    bilinearForm.getLocalMatrix(bilinearFormMatrix);

    // compute Bu - f (we do f-= Bu so the output is in rhsElement)
    bilinearFormMatrix.mmv(solutionElement,rhsElement);

    // Solve for the Riesz lift and then compute the residual
    BlockVector<FieldVector<double,1> > tmpVector(rhsElement);
    Cholesky<Matrix<FieldMatrix<double,1,1>>> cholesky(innerProductMatrix);
    cholesky.apply(tmpVector);

    // compute the scalar product of the following two vectors
    return tmpVector * rhsElement;
  }

/**
 * \brief Computation of a posteriori error in (u,theta)
 *
 * \param bilinearForm        the bilinear form
 * \param innerProduct        the inner product
 * \param solution           the computed solution
 * \param rhs                 the right-hand side
 */
  template <class BilinearForm,class InnerProduct,
            class VectorType>
  double ErrorTools::aPosterioriError(BilinearForm& bilinearForm,
                                      InnerProduct& innerProduct,
                                      const VectorType& solution,
                                      const VectorType& rhs)
  {
    using namespace boost::fusion;
    using namespace Dune::detail;

    typedef typename std::tuple_element<0,typename BilinearForm::SolutionSpaces>::type::GridView GridView;
    const GridView gridView = std::get<0>(bilinearForm.getSolutionSpaces()).gridView();

    typedef typename BilinearForm::SolutionSpaces SolutionSpaces;
    typedef typename BilinearForm::TestSpaces EnrichedTestspaces;

    typedef typename result_of::as_vector<typename result_of::transform<SolutionSpaces,getLocalView>::type>::type SolutionLocalViews;
    typedef typename result_of::as_vector<typename result_of::transform<EnrichedTestspaces,getLocalView>::type>::type TestLocalViews;

    SolutionLocalViews localViewsSolution = as_vector(transform(bilinearForm.getSolutionSpaces(),getLocalView()));
    TestLocalViews localViewsTest = as_vector(transform(bilinearForm.getTestSpaces(),getLocalView()));

    // We get the local index sets of the test spaces
    auto testLocalIndexSets = as_vector(transform(bilinearForm.getTestSpaces(),
                                                  getLocalIndexSet()));
    // We get the local index sets of the solution spaces
    auto solutionLocalIndexSets
            = as_vector(transform(bilinearForm.getSolutionSpaces(),
                                  getLocalIndexSet()));

    // Variable where we compute the residual
    double res = 0.;

    for(const auto& e : elements(gridView))
    {
      // Bind localViews and localIndexSets
      for_each(localViewsTest, applyBind<decltype(e)>(e));
      for_each(localViewsSolution, applyBind<decltype(e)>(e));
      for_each(zip(testLocalIndexSets, localViewsTest),
               make_fused_procedure(bindLocalIndexSet()));
      for_each(zip(solutionLocalIndexSets, localViewsSolution),
               make_fused_procedure(bindLocalIndexSet()));

      res += aPosterioriErrorSquareElement(bilinearForm,
                                           innerProduct,
                                           localViewsTest,
                                           localViewsSolution,
                                           testLocalIndexSets,
                                           solutionLocalIndexSets,
                                           solution,
                                           rhs);
   }

   return std::sqrt(res);

  }

/**
 * \brief Mark elements for refinement according to Dörfler's strategy
 *
 * This mean, given a ratio $\theta \in (0,1]$, the set
 * $\mathcal M \subset \Omega_h$ of marked elements satisfies
 * $$\operatorname{err}(u_h, \mathcal M) \geq \operatorname{err}(u_h, \Omega_h).$$
 *
 * \param ratio  the marking ratio $\theta \in (0,1]$
 * \param bilinearForm
 * \param innerProduct
 * \param solution  FE solution
 * \param rhs  Rhs of the problem in enriched test space
 *
 * \return estimate for global a posteriori error
 */
  template <class Grid, class BilinearForm, class InnerProduct,
            class APosterioriInnerProduct, class LinearForm,
            class RhsFunction, class VectorType>
  double ErrorTools::DoerflerMarking(
      Grid& grid,
      double ratio,
      BilinearForm& bilinearForm,
      InnerProduct& innerProduct,
      APosterioriInnerProduct& aPosterioriInnerProduct,
      LinearForm& linearForm,
      const RhsFunction& f,
      const VectorType& solution,
      const VectorType& rhs,
      double splitRatio
                                    )
  {
    using namespace boost::fusion;
    using namespace Dune::detail;

    assert(splitRatio >= 0 && splitRatio <= 1);
    using GridView
        = typename std::tuple_element<0,typename BilinearForm::SolutionSpaces>
                        ::type::GridView;
    static_assert(std::is_same<typename GridView::Grid, Grid>::value,
        "Type mismatch between Grid and Grid of BilinearForm!");
    using Entity = typename GridView::template Codim<0>::Entity;
    using EntitySeed = typename Entity::EntitySeed;

    const GridView gridView
        = std::get<0>(bilinearForm.getSolutionSpaces()).gridView();

    using SolutionSpaces = typename BilinearForm::SolutionSpaces;
    using EnrichedTestspaces = typename BilinearForm::TestSpaces;

    using SolutionLocalViews
        = typename result_of::as_vector<typename
            result_of::transform<SolutionSpaces,getLocalView>::type>::type;
    using TestLocalViews
        = typename result_of::as_vector<typename
            result_of::transform<EnrichedTestspaces,getLocalView>::type>::type;

    SolutionLocalViews localViewsSolution
        = as_vector(transform(bilinearForm.getSolutionSpaces(), getLocalView()));
    TestLocalViews localViewsTest
        = as_vector(transform(bilinearForm.getTestSpaces(), getLocalView()));

    // We get the local index sets of the test spaces
    auto testLocalIndexSets = as_vector(transform(bilinearForm.getTestSpaces(),
                                                  getLocalIndexSet()));
    // We get the local index sets of the solution spaces
    auto solutionLocalIndexSets
            = as_vector(transform(bilinearForm.getSolutionSpaces(),
                                  getLocalIndexSet()));

    std::vector<std::tuple<EntitySeed, double>> errorEstimates;
    errorEstimates.reserve(gridView.size(0));
    for(const auto& e : elements(gridView))
    {
      // Bind localViews and localIndexSets
      for_each(localViewsTest, applyBind<decltype(e)>(e));
      for_each(localViewsSolution, applyBind<decltype(e)>(e));
      for_each(zip(testLocalIndexSets, localViewsTest),
               make_fused_procedure(bindLocalIndexSet()));
      for_each(zip(solutionLocalIndexSets, localViewsSolution),
               make_fused_procedure(bindLocalIndexSet()));

      // Now we compute the error inside the element
      double elementError = 0;
      if (splitRatio > 0)
      {
        elementError += splitRatio*aPosterioriErrorSquareElement(
                                bilinearForm,
                                innerProduct,
                                localViewsTest,
                                localViewsSolution,
                                testLocalIndexSets,
                                solutionLocalIndexSets,
                                solution,
                                rhs);
      }
      if (splitRatio < 1)
      {
        elementError += (1-splitRatio)*aPosterioriL2ErrorSquareElement(
                                aPosterioriInnerProduct,
                                linearForm,
                                f,
                                localViewsSolution,
                                solutionLocalIndexSets,
                                solution);
      }
      errorEstimates.emplace_back(e.seed(), elementError);
    }
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

    return std::sqrt(errorSquared);
  }

} // end namespace Dune

#endif // DUNE_DPG_ERROR_TOOLS
