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
    template <class LocalView,class VolumeTerms>
    double computeL2errorElement(const LocalView& ,
                                 BlockVector<FieldVector<double,1> >& ,
                                 VolumeTerms&& );

    template <class FEBasis,class VolumeTerms>
    double computeL2error(const FEBasis& ,
                          BlockVector<FieldVector<double,1> >& ,
                          VolumeTerms&&);

    template<class GridType> void hRefinement(GridType& grid);

    template <class BilinearForm,class InnerProduct,class VectorType>
    double aPosterioriError(BilinearForm& ,
                            InnerProduct& ,
                            VectorType& ,
                            VectorType& ,
                            VectorType& );

  };

//*******************************************************************

/**
 * \brief Returns the computation in a given element of the L2 error
          between the exact solution uRef and the fem solution u.
 *
 * \param localView      the local view of the element
 * \param u              the vector containing the computed solution
 * \param uRef           the expression for the exact solution.
 */
  template <class LocalView,class VolumeTerms>
  double ErrorTools::computeL2errorElement(const LocalView& localView,
                                           BlockVector<FieldVector<double,1> >& u,
                                           VolumeTerms&& uRef)
  {

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    const int dim = Element::dimension;
    auto geometry = element.geometry();

    // Get set of shape functions for this element
    const auto& localFiniteElement = localView.tree().finiteElement();

    const unsigned int quadratureOrder = 10; // Remark: the quadrature order has to be an even number
    const QuadratureRule<double, dim>& quad =
        QuadratureRules<double, dim>::rule(element.type(), quadratureOrder);

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
      const FieldVector<double,dim>& mapQuadPos = geometry.global(quadPos) ;  // we get the global coordinate of quadPos

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
 * \param feBasis        the finite element basis
 * \param u              the vector containing the computed solution
 * \param uRef           the expression for the exact solution.
 */
  template <class FEBasis,class VolumeTerms>
  double ErrorTools::computeL2error(const FEBasis& feBasis,
                                    BlockVector<FieldVector<double,1> >& u,
                                    VolumeTerms&& uRef)
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

    // Position of the coefficients
    int posBegin = 0;
    int posEnd   = 0;

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
      errorElement_[indexElement] = computeL2errorElement(localView,uElement,uRef);
      errSquare += errorElement_[indexElement]*errorElement_[indexElement];
    }

    return std::sqrt(errSquare);
  }

/**
 * \brief Computation of a posteriori error in (u,theta)
 *
 * \param bilinearForm        the bilinear form
 * \param innerProduct        the inner product
 * \param uSolution           the computed solution u
 * \param thetaSolution       the computed solution theta
 * \param rhs                 the right-hand side
 */
  template <class BilinearForm,class InnerProduct,class VectorType>
  double ErrorTools::aPosterioriError(BilinearForm& bilinearForm,
                                      InnerProduct& innerProduct,
                                      VectorType&   uSolution,
                                      VectorType&   thetaSolution,
                                      VectorType&   rhs)
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
      for_each(localViewsTest, applyBind<decltype(e)>(e));
      for_each(localViewsSolution, applyBind<decltype(e)>(e));
      for_each(zip(testLocalIndexSets, localViewsTest),
               make_fused_procedure(bindLocalIndexSet()));
      for_each(zip(solutionLocalIndexSets, localViewsSolution),
               make_fused_procedure(bindLocalIndexSet()));

      // We take the coefficients of our solution (u,theta) that correspond to the current element e. They are stored in uLocal and thetaLocal.
      // dof of the finite element inside the element
      size_t dofElementU = boost::fusion::at_c<0>(localViewsSolution).size();
      size_t dofElementTheta = boost::fusion::at_c<1>(localViewsSolution).size();

      // We extract the solution on the current element
      BlockVector<FieldVector<double,1> > solutionElement(dofElementU+dofElementTheta);
      // We take the coefficients of uSolution that correspond to the current element e. They are stored in the second part of solutionElement.
      for (size_t i=0; i<dofElementU; i++)
      {
        solutionElement[i] = uSolution[ at_c<0>(solutionLocalIndexSets).index(i)[0] ];
      }
      // We take the coefficients of thetaSolution that correspond to the current element e. They are stored in the first part of solutionElement.
      for (size_t i=0; i<dofElementTheta; i++)
      {
        solutionElement[i+dofElementU] = thetaSolution[ at_c<1>(solutionLocalIndexSets).index(i)[0] ];
      }

      // We get rhsLocal
      // dof of the finite element test space inside the element
      size_t dofElementTest = boost::fusion::at_c<0>(localViewsTest).size();
      // We extract the rhs of the current element
      BlockVector<FieldVector<double,1> > rhsElement(dofElementTest);
      for (size_t i=0; i<dofElementTest; i++)
      {
        rhsElement[i] = rhs[ at_c<0>(testLocalIndexSets).index(i)[0] ];
      }

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
