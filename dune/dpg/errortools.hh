#ifndef DUNE_DPG_ERROR_TOOLS
#define DUNE_DPG_ERROR_TOOLS


#include <iostream>

#include <vector>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dune/dpg/cholesky.hh>


namespace Dune {


  //*******************************************************************
  class ErrorTools
  {
    std::vector<double> errorElement_;

  public:
    ErrorTools() {};
    template <unsigned int subsamples, class LocalView, class VolumeTerms>
    double computeL2errorElement(const LocalView& ,
                                 BlockVector<FieldVector<double,1> >& ,
                                 VolumeTerms&&,
                                 unsigned int = 5);

    template <unsigned int subsamples, class FEBasis, class VolumeTerms>
    double computeL2error(const FEBasis& ,
                          BlockVector<FieldVector<double,1> >& ,
                          VolumeTerms&&,
                          unsigned int = 5);

    template <class BilinearForm,class InnerProduct,class VectorType>
    double aPosterioriErrorSquareElement(BilinearForm& ,
                                      InnerProduct& ,
                                      VectorType& ,
                                      VectorType& );

    template <class BilinearForm, class InnerProduct, class VectorType>
    double aPosterioriError(BilinearForm& ,
                            InnerProduct& ,
                            const VectorType& ,
                            const VectorType& );

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
  double ErrorTools::computeL2errorElement(const LocalView& localView,
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

    return std::sqrt(errSquare);
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

     //   Number of elements
    const int codimElement = 0;   //Element: codim = 0
    int nElement = gridView.indexSet().size(codimElement);

    // we update the size of errorElement_ where we will store the errors per element
    errorElement_.resize(nElement);

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
      errorElement_[indexElement] = computeL2errorElement<subsamples>
                                      (localView,uElement,uRef,quadratureOrder);
      errSquare += errorElement_[indexElement]*errorElement_[indexElement];
    }

    return std::sqrt(errSquare);
  }



  /**
 * \brief Computation of a posteriori error in (u,theta)
 *
 * \param bilinearForm        the bilinear form
 * \param innerProduct        the inner product
 * \param solution           the computed solution
 * \param rhs                 the right-hand side
 */
  template <class BilinearForm,class InnerProduct,class VectorType>
  double ErrorTools::aPosterioriErrorSquareElement(BilinearForm& bilinearForm,
                                      InnerProduct& innerProduct,
                                      VectorType& solution,
                                      VectorType& rhs)
  {

  }

  /**
 * \brief Computation of a posteriori error in (u,theta)
 *
 * \param bilinearForm        the bilinear form
 * \param innerProduct        the inner product
 * \param solution           the computed solution
 * \param rhs                 the right-hand side
 */
  template <class BilinearForm,class InnerProduct,class VectorType>
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
      res += tmpVector * rhsElement;

   }

   return std::sqrt(res);

  }

} // end namespace Dune

#endif // DUNE_DPG_ERROR_TOOLS
