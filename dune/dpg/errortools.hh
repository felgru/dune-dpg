#ifndef DUNE_DPG_ERROR_TOOLS
#define DUNE_DPG_ERROR_TOOLS


#include <iostream>

#include <vector>

#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dune/istl/io.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/grid/uggrid.hh>   // for triangular meshes that are locally adaptive
#include <dune/grid/utility/structuredgridfactory.hh> // for triangular meshes that are locally adaptive


namespace Dune {


  //*******************************************************************
  class ErrorTools
  {
    double errorTol_;
    std::vector<double> errorElement_ ;

  public:
    ErrorTools(double);
    template <class LocalView,class VolumeTerms>
    double computeL2errorElement(
                                 const LocalView& ,
                                 BlockVector<FieldVector<double,1> >& ,
                                 VolumeTerms&
                                );

    template <class FEBasis,class VolumeTerms>
    double computeL2error(
                          const FEBasis& ,
                          BlockVector<FieldVector<double,1> >& ,
                          VolumeTerms&
                         );

    template<class GridType> void hRefinement( GridType& grid );

    template <class BilinearForm,class InnerProduct,class VectorType>
    double aPosterioriError(
                            BilinearForm& ,
                            InnerProduct& ,
                            VectorType& ,
                            VectorType& ,
                            VectorType&
                           );

  };

//*******************************************************************
  ErrorTools::ErrorTools (double tol)
  {
    errorTol_ = tol;
  }

//*******************************************************************

  // computatiton of an L2 error in an element
  template <class LocalView,class VolumeTerms>
  double ErrorTools::computeL2errorElement(
                                const LocalView& localView,
                                BlockVector<FieldVector<double,1> >& u,
                                VolumeTerms& uRef
                              )
  {

    // Get the grid element from the local FE basis view
    typedef typename LocalView::Element Element;
    const Element& element = localView.element();

    const int dim = Element::dimension;
    auto geometry = element.geometry();

    // Get set of shape functions for this element
    const auto& localFiniteElement = localView.tree().finiteElement();

    // Get a quadrature rule
    int order = 10; // Remark: the quadrature order has to be an even number (nombre pair)!
    const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element.type(), order);

    // Variables employed in the loop
    double errSquare = 0 ;     // we store here the square of the error
    double uQuad = 0 ;         // we store here the value of u at a quadrature point

    // Loop over all quadrature points
    for (size_t pt=0; pt < quad.size(); pt++) {

      //Restart value uQuad
      uQuad = 0 ;

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();
      // Position of the current quadrature point in the current element
      const FieldVector<double,dim>& mapQuadPos = geometry.global(quadPos) ;  // we get the global coordinate of quadPos

      // The multiplicative factor in the integral transformation formula
      const double integrationElement = element.geometry().integrationElement(quadPos);

      // Evaluate all shape function values at mapQuadPos (we do this through quadPos, which is the point correpondind to quadPos in the reference element)
      std::vector<FieldVector<double,1> > shapeFunctionValues;
      localFiniteElement.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

      // Evaluation of u at the point mapQuadPos, which is the point correpondind to quadPos in the reference element
      for(int i=0; i<shapeFunctionValues.size(); i++)
      {
        uQuad += shapeFunctionValues[i]*u[i] ;
      }

      // Value of uRef at mapQuadPos
      double uExactQuad = std::get<0>(uRef)(mapQuadPos);

      // we add the error at the quadrature point
      errSquare += (uQuad - uExactQuad)*(uQuad - uExactQuad) * quad[pt].weight() * integrationElement ;

    }

    return std::sqrt(errSquare) ;
  }

//*******************************************************************
// computation in the whole mesh of the L2 error
// between the exact solution uRef and the fem solution u
  template <class FEBasis,class VolumeTerms>
  double ErrorTools::computeL2error(
                        const FEBasis& feBasis,
                        BlockVector<FieldVector<double,1> >& u,
                        VolumeTerms& uRef
                        )
  {

    // Get the grid view from the finite element basis
    typedef typename FEBasis::GridView GridView;
    GridView gridView = feBasis.gridView();

     //   Number of elements
    const int codimElement = 0 ;   //Element: codim = 0
    int nElement = gridView.indexSet().size(codimElement) ;

    // we update the size of errorElement_ where we will store the errors per element
    errorElement_.resize(nElement);

    // Variables where we will store the errors
    double errSquare = 0.;

    // Get index set of the finite element basis
    auto basisIndexSet = feBasis.indexSet();
    // A view on the FE basis on a single element
    auto localView = feBasis.localView();
    auto localIndexSet = basisIndexSet.localIndexSet();

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
      //Remark syntax: The global index of the i-th vertex (because it is codim [0]) of the element 'e'
      //auto row = localIndexSet.index(i)[0]; //equivalent to indexSet.subIndex(∗it, j , dim);

      // Now we compute the error inside the element
      errorElement_[indexElement] = computeL2errorElement(localView,uElement,uRef);
      errSquare += errorElement_[indexElement]*errorElement_[indexElement] ;
    }

    return std::sqrt(errSquare);
  }

//*******************************************************************
  template<class GridType> void ErrorTools::hRefinement( GridType& grid )
  {

  //typedef GridType::LeafGridView GridView;
  auto gridView = grid->leafGridView();
  // GridView::Codim<0>::Iterator it = gridView.begin<0>();
  // GridView::Codim<0>::Iterator endIt = gridView.end<0>();
  //for (; it != endIt; ++it)
  for(const auto& e : elements(gridView))
  {
    //int idx = mapper.index(*it); // get element index  // Remark, the old method for this is mapper.map(*it);

    int indexElement = gridView.indexSet().index(e);
    std::cout << "element with index " << indexElement << std::endl ;
    if ( errorElement_[indexElement] < errorTol_)
    {
      grid->mark(0, e); // index 0: we do not mark for refinement nor coarsening
      std::cout << "Error element = " << errorElement_[indexElement] << " --> not h-refinement " << std::endl ;
    }
    else if ( errorElement_[indexElement] >= errorTol_)
    {
      grid->mark(1, e); // index 1: mark for refinement
      std::cout << "Error element = " << errorElement_[indexElement] << " -->  h-refinement " << std::endl ;
    }

  }

  grid->preAdapt();
  grid->adapt();
  grid->postAdapt();

  }

//*******************************************************************
  // computation of a posteriori error in (u,theta)
  template <class BilinearForm,class InnerProduct,class VectorType>
  double ErrorTools::aPosterioriError(
                        BilinearForm& bilinearForm,
                        InnerProduct& innerProduct,
                        VectorType& uSolution,
                        VectorType& thetaSolution,
                        VectorType& rhs
                        )
  {
    using namespace boost::fusion;
    using namespace Dune::detail;

    typedef typename std::tuple_element<0,typename BilinearForm::SolutionSpaces>::type::GridView GridView;
    const GridView gridView = std::get<0>(bilinearForm.getSolutionSpaces()).gridView() ;

    typedef typename BilinearForm::SolutionSpaces SolutionSpaces;
    typedef typename BilinearForm::TestSpaces EnrichedTestspaces;

    typedef typename result_of::as_vector<typename result_of::transform<SolutionSpaces,getLocalView>::type>::type SolutionLocalView;
    typedef typename result_of::as_vector<typename result_of::transform<EnrichedTestspaces,getLocalView>::type>::type TestLocalView;

    SolutionLocalView localViewSolution = as_vector(transform(bilinearForm.getSolutionSpaces(),getLocalView()));
    TestLocalView localViewTest = as_vector(transform(bilinearForm.getTestSpaces(),getLocalView()));

    //We get the index set of the test spaces: variable testLocalIndexSet
    auto testBasisIndexSet = as_vector(transform(bilinearForm.getTestSpaces(), getIndexSet()));
    auto testLocalIndexSet = as_vector(transform(testBasisIndexSet,getLocalIndexSet()));
    for_each(zip(testLocalIndexSet, localViewTest),make_fused_procedure(bindLocalIndexSet()));
    //We get the index set of the solution spaces: variable solutionLocalIndexSet
    auto solutionBasisIndexSet = as_vector(transform(bilinearForm.getSolutionSpaces(), getIndexSet()));
    auto solutionLocalIndexSet = as_vector(transform(solutionBasisIndexSet,getLocalIndexSet()));
    for_each(zip(solutionLocalIndexSet, localViewSolution),make_fused_procedure(bindLocalIndexSet()));

    // Variable where we compute the residual
    double res = 0.;

    for(const auto& e : elements(gridView))
    {
      for_each(localViewTest, applyBind<decltype(e)>(e));
      for_each(localViewSolution, applyBind<decltype(e)>(e));

      // We take the coefficients of our solution (u,theta) that correspond to the current element e. They are stored in uLocal anf thetaLocal.
      // dof of the finite element inside the element
      size_t dofElementTheta = boost::fusion::at_c<0>(localViewSolution)->size();
      size_t dofElementU = boost::fusion::at_c<1>(localViewSolution)->size();

      // We extract the solution on the current element
      BlockVector<FieldVector<double,1> > solutionElement(dofElementTheta+dofElementU);
      // We take the coefficients of thetaSolution that correspond to the current element e. They are stored in the first part of solutionElement.
      for (size_t i=0; i<dofElementTheta; i++)
      {
        solutionElement[i] = thetaSolution[ at_c<0>(solutionLocalIndexSet)->index(i)[0] ];
      }
      // We take the coefficients of uSolution that correspond to the current element e. They are stored in the second part of solutionElement.
      for (size_t i=0; i<dofElementU; i++)
      {
        solutionElement[i+dofElementTheta] = uSolution[ at_c<1>(solutionLocalIndexSet)->index(i)[0] ];
      }

      // We get rhsLocal
      // dof of the finite element test space inside the element
      size_t dofElementTest = boost::fusion::at_c<0>(localViewTest)->size();
      // We extract the rhs of the current element
      BlockVector<FieldVector<double,1> > rhsElement(dofElementTest);
      for (size_t i=0; i<dofElementTest; i++)
      {
        rhsElement[i] = rhs[ at_c<0>(testLocalIndexSet)->index(i)[0] ];
      }

      // We grab the inner product matrix in the innerProductMatrix variable (IP)
      Matrix<FieldMatrix<double,1,1> > innerProductMatrix;
      innerProduct.bind(localViewTest);
      innerProduct.getLocalMatrix(innerProductMatrix);

      // We grab the bilinear form matrix in the bilinearFormMatirx variable
      Matrix<FieldMatrix<double,1,1> > bilinearFormMatirx;
      bilinearForm.bind(localViewTest, localViewSolution);
      bilinearForm.getLocalMatrix(bilinearFormMatirx);

      // compute Bu - f (we do f-= Bu so the output is in rhsElement)
      bilinearFormMatirx.mmv(solutionElement,rhsElement);

      //Solve the system
      // We first copy the matrix because it is going to be destroyed by the Cholesky
      Matrix<FieldMatrix<double,1,1> > innerProductMatrixCholesky = innerProductMatrix;

      Functions::Cholesky<Matrix<FieldMatrix<double,1,1> > > cholesky(innerProductMatrixCholesky);
      cholesky.apply(rhsElement); // the solution is overriten in rhsElement

      // computation of the residual
      BlockVector<FieldVector<double,1> > tmpVector(rhsElement.size());
      innerProductMatrix.mv(rhsElement,tmpVector);
      res += rhsElement * tmpVector; // remark: the operation * in vectors of Dune computes the scalar product of two vectors

   }

   return std::sqrt(res);

  }

} // end namespace Dune

#endif // DUNE_DPG_ERROR_TOOLS
