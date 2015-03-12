#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <vector>

#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/common/intersection.hh> //TODO necessary?

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/istl/umfpack.hh>

#define DUNE_FINAL final

#include <dune/functions/functionspacebases/q2tracenodalbasis.hh>
#include <dune/functions/functionspacebases/qkdiscontinuousnodalbasis.hh>

#include <dune/dpg/functionhelper.hh>

using namespace Dune;






// Compute the stiffness matrix for a single element
template <class LocalViewTrace, class LocalViewInterior, class LocalViewTest, class GridView, class MatrixType>
void getLocalMatrix(const LocalViewTrace& localViewTrace,
                    const LocalViewInterior& localViewInterior,
                    const LocalViewTest& localViewTest,
                    const GridView& gridView,
                    MatrixType& elementMatrix)
{
  FieldVector<double, 2> beta={1,1};
  double c=0;

  // Get the grid element from the local FE basis view
  typedef typename LocalViewTest::Element Element;
  const Element& element = localViewTest.element();

  const int dim = Element::dimension;
  auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& localFiniteElementTest = localViewTest.tree().finiteElement();
  const auto& localFiniteElementTrace = localViewTrace.tree().finiteElement();
  const auto& localFiniteElementInterior = localViewInterior.tree().finiteElement();

  const int nTest(localFiniteElementTest.localBasis().size());
  const int nTrace(localFiniteElementTrace.localBasis().size());
  const int nInterior(localFiniteElementInterior.localBasis().size());

  // Set all matrix entries to zero
  elementMatrix.setSize(nTest, nTest+nTrace+nInterior);
  elementMatrix = 0;      // fills the entire matrix with zeroes

  ////////////////////////////
  // Assemble interior terms
  ////////////////////////////

  // Get a quadrature rule
  //int order = 2*(dim-1); //TODO!!!!!!!!!
  int order = 20;
  const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (size_t pt=0; pt < quad.size(); pt++) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The transposed inverse Jacobian of the map from the reference element to the element
    const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPos);

    ////////////////////////////
    // Test Functions:
    ////////////////////////////
    // values of the shape functions
    std::vector<FieldVector<double,1> > valuesTest;
    localFiniteElementTest.localBasis().evaluateFunction(quadPos, valuesTest);

    // The gradients of the shape functions on the reference element
    std::vector<FieldMatrix<double,1,dim> > referenceGradientsTest;
    localFiniteElementTest.localBasis().evaluateJacobian(quadPos, referenceGradientsTest);

    // Compute the shape function gradients on the real element
    std::vector<FieldVector<double,dim> > gradientsTest(referenceGradientsTest.size());
    for (size_t i=0; i<gradientsTest.size(); i++)
      jacobian.mv(referenceGradientsTest[i][0], gradientsTest[i]);

    ////////////////////////////
    // Interior Trial Functions
    ////////////////////////////
    // values of the shape functions
    std::vector<FieldVector<double,1> > valuesInterior;
    localFiniteElementInterior.localBasis().evaluateFunction(quadPos, valuesInterior);

    // Compute the actual matrix entries
    for (size_t i=0; i<nTest; i++)
    {
      for (size_t j=0; j<nTest; j++ )
      {
        elementMatrix[i][j] += (valuesTest[i] * valuesTest[j]) * quad[pt].weight() * integrationElement;
        elementMatrix[i][j] += (beta*gradientsTest[i]) * (beta*gradientsTest[j]) * quad[pt].weight() * integrationElement;
      }
      for (size_t j=0; j<nInterior; j++)
      {
        elementMatrix[i][j+nTest+nTrace] += (valuesTest[i] * valuesInterior[j])* c * quad[pt].weight() * integrationElement;
        elementMatrix[i][j+nTest+nTrace] += (beta*gradientsTest[i])*(-1.0) * valuesInterior[j] * quad[pt].weight() * integrationElement;
      }
    }
  }


  ////////////////////////////
  // Assemble boundary terms
  ////////////////////////////

  for (auto&& intersection : intersections(gridView, element))
  {
    const QuadratureRule<double, dim-1>& quadFace = QuadratureRules<double, dim-1>::rule(intersection.type(), order);
    // Loop over all quadrature points
    for (size_t pt=0; pt < quadFace.size(); pt++) {

    // Position of the current quadrature point in the reference element (face!)
    const FieldVector<double,dim-1>& quadFacePos = quadFace[pt].position();

    // The multiplicative factor in the integral transformation formula moltiplied with outer normal
    const FieldVector<double,dim>& integrationOuterNormal = intersection.integrationOuterNormal(quadFacePos);

                // position of the quadrature point within the element
    const FieldVector<double,dim> elementQuadPos = intersection.geometryInInside().global(quadFacePos);


    ////////////////////////////
    // Test Functions:
    ////////////////////////////
    // values of the shape functions
    std::vector<FieldVector<double,1> > valuesTest;
    localFiniteElementTest.localBasis().evaluateFunction(elementQuadPos, valuesTest);

    ////////////////////////////
    // Trace Functions
    ////////////////////////////
    // values of the shape functions
    std::vector<FieldVector<double,1> > valuesTrace;
    localFiniteElementTrace.localBasis().evaluateFunction(elementQuadPos, valuesTrace);

    // Compute the actual matrix entries
    for (size_t i=0; i<nTest; i++)
    {
      for (size_t j=0; j<nTrace; j++ )
      {
        elementMatrix[i][j+nTest] += ((beta*integrationOuterNormal) * valuesTest[i] * valuesTrace[j]) * quadFace[pt].weight();
      }
    }
    }
  }

}



// Compute the source term for a single element
template <class LocalViewTest>
void getVolumeTerm( const LocalViewTest& localViewTest,
                    BlockVector<FieldVector<double,1> >& localRhs,
                    const Dune::VirtualFunction<FieldVector<double,LocalViewTest::Element::dimension>, double>* volumeTerm)
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

  // A quadrature rule
  int order = 20; //TODO!!!!!!
  const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element.type(), order);


  // Loop over all quadrature points
  for ( size_t pt=0; pt < quad.size(); pt++ ) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = element.geometry().integrationElement(quadPos);

    double functionValue;
    volumeTerm->evaluate(element.geometry().global(quadPos), functionValue);

    // Evaluate all shape function values at this point
    std::vector<FieldVector<double,1> > shapeFunctionValues;
    localFiniteElementTest.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

    // Actually compute the vector entries
    for (size_t i=0; i<localRhs.size(); i++)
      localRhs[i] += shapeFunctionValues[i] * functionValue * quad[pt].weight() * integrationElement;

  }

}

// Get the occupation pattern of the stiffness matrix TODO!!!!!!!!
template <class FEBasisTrace, class FEBasisInterior, class FEBasisTest>
void getOccupationPattern(const FEBasisTrace& feBasisTrace,
                            const FEBasisInterior& feBasisInterior,
                            const FEBasisTest& feBasisTest,
                            MatrixIndexSet& nb)
{
  // Total number of degrees of freedom
  auto n = feBasisTrace.subIndexCount()+feBasisInterior.subIndexCount()+feBasisTest.subIndexCount();

  nb.resize(n, n);

  // A view on the FE basis on a single element
  typename FEBasisTest::LocalView localViewTest(&feBasisTest);
  typename FEBasisTrace::LocalView localViewTrace(&feBasisTrace);
  typename FEBasisInterior::LocalView localViewInterior(&feBasisInterior);

  // Loop over all leaf elements
  auto it    = feBasisTest.gridView().template begin<0>();
  auto endIt = feBasisTest.gridView().template end<0>  ();
  auto itTrace    = feBasisTrace.gridView().template begin<0>();
  auto itInterior    = feBasisInterior.gridView().template begin<0>();

  int offsetTrace(feBasisTest.subIndexCount());
  int offsetInterior(offsetTrace+feBasisTrace.subIndexCount());

  for (; it!=endIt; ++it, ++itTrace, ++itInterior)
  {
    // Bind the local FE basis view to the current element
    localViewTest.bind(*it);
    localViewTrace.bind(*itTrace);
    localViewInterior.bind(*itInterior);

        // There is a matrix entry a_ij if the i-th and j-th vertex are connected in the grid
        for (size_t i=0; i<localViewTest.tree().localSize(); i++) {

            for (size_t j=0; j<localViewTest.tree().localSize(); j++) {

                auto iIdx = localViewTest.tree().globalIndex(i)[0];
                auto jIdx = localViewTest.tree().globalIndex(j)[0];

                // Add a nonzero entry to the matrix
                nb.add(iIdx, jIdx);

            }
            for (size_t j=0; j<localViewTrace.tree().localSize(); j++) {

                auto iIdx = localViewTest.tree().globalIndex(i)[0];
                auto jIdx = localViewTrace.tree().globalIndex(j)[0];

                // Add a nonzero entry to the matrix
                nb.add(iIdx, jIdx+offsetTrace);
                nb.add(jIdx+offsetTrace, iIdx);

            }
            for (size_t j=0; j<localViewInterior.tree().localSize(); j++) {

                auto iIdx = localViewTest.tree().globalIndex(i)[0];
                auto jIdx = localViewInterior.tree().globalIndex(j)[0];

                // Add a nonzero entry to the matrix
                nb.add(iIdx, jIdx+offsetInterior);
                nb.add(jIdx+offsetInterior, iIdx);
            }
        }
    }
    for (size_t i=0; i<n; i++)
    {
      nb.add(i, i);
    }

}


/** \brief Assemble the Laplace stiffness matrix on the given grid view */
template <class FEBasisTrace, class FEBasisInterior, class FEBasisTest>
void assembleSystem(const FEBasisTrace& feBasisTrace,
                    const FEBasisInterior& feBasisInterior,
                    const FEBasisTest& feBasisTest,
                    BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                    BlockVector<FieldVector<double,1> >& rhs,
                    const VirtualFunction<Dune::FieldVector<double,FEBasisTest::GridView::dimension>,double>* volumeTerm)
{
  // Get the grid view from the finite element basis
  typedef typename FEBasisTest::GridView GridView;
  GridView gridView = feBasisTest.gridView();

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  MatrixIndexSet occupationPattern;
  getOccupationPattern<FEBasisTrace, FEBasisInterior, FEBasisTest>(feBasisTrace, feBasisInterior, feBasisTest, occupationPattern);

  // ... and give it the occupation pattern we want.
  occupationPattern.exportIdx(matrix);

  // set rhs to correct length -- the total number of basis vectors in the basis
  rhs.resize(feBasisTrace.subIndexCount()+feBasisInterior.subIndexCount()+feBasisTest.subIndexCount());

  // Set all entries to zero
  matrix = 0;
  rhs = 0;

  // A view on the FE basis on a single element
  typename FEBasisTrace::LocalView localViewTrace(&feBasisTrace);
  typename FEBasisInterior::LocalView localViewInterior(&feBasisInterior);
  typename FEBasisTest::LocalView localViewTest(&feBasisTest);

  // A loop over all elements of the grid
  auto it    = gridView.template begin<0>();
  auto endIt = gridView.template end<0>  ();

  for( ; it != endIt; ++it ) {

    // Bind the local FE basis view to the current element
    localViewTrace.bind(*it);
    localViewInterior.bind(*it);
    localViewTest.bind(*it);

    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    getLocalMatrix(localViewTrace, localViewInterior, localViewTest, gridView, elementMatrix);

    const int nTest(localViewTest.tree().finiteElement().localBasis().size());
    const int nTrace(localViewTrace.tree().finiteElement().localBasis().size());
    const int nInterior(localViewInterior.tree().finiteElement().localBasis().size());

    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i=0; i<nTest; i++)
    {
      // The global index of the i-th local degree of freedom of the element 'it'
      auto row = localViewTest.tree().globalIndex(i)[0];

      for (size_t j=0; j<nTest; j++ )
      {
        // The global index of the j-th local degree of freedom of the element 'it'
        auto col = localViewTest.tree().globalIndex(j)[0];
        matrix[row][col] += elementMatrix[i][j];
      }
      for (size_t j=0; j<nTrace; j++ )
      {
        // The global index of the j-th local degree of freedom of the element 'it'
        auto col = localViewTrace.tree().globalIndex(j)[0]
                    +feBasisTest.subIndexCount();
        matrix[row][col] += elementMatrix[i][j+nTest];
        matrix[col][row] += elementMatrix[i][j+nTest];
      }
      for (size_t j=0; j<nInterior; j++ )
      {
        // The global index of the j-th local degree of freedom of the element 'it'
        auto col = localViewInterior.tree().globalIndex(j)[0]
                    +feBasisTest.subIndexCount()
                    +feBasisTrace.subIndexCount();
        matrix[row][col] += elementMatrix[i][j+nTest+nTrace];
        matrix[col][row] += elementMatrix[i][j+nTest+nTrace];
      }
    }

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> > localRhs;
    getVolumeTerm(localViewTest, localRhs, volumeTerm);

     for (size_t i=0; i<localRhs.size(); i++) {

      // The global index of the i-th vertex of the element 'it'
      auto row = localViewTest.tree().globalIndex(i)[0];
      rhs[row] += localRhs[i];

    }

  }

}

// The identity function in R^dim.
// This function is used to find the positions of the Lagrange nodes of the FE space
template <int dim>
struct Identity
    : public VirtualFunction<FieldVector<double,dim>, FieldVector<double,dim> >
{
    void evaluate(const FieldVector<double,dim>& in, FieldVector<double,dim>& out) const {
        out = in;
    }
};


// TODO !!!!!!!!!!
// This method marks all vertices on the boundary of the grid.
// In our problem these are precisely the Dirichlet nodes.
// The result can be found in the 'dirichletNodes' variable.  There, a bit
// is set precisely when the corresponding vertex is on the grid boundary.
template <class FEBasis>
void boundaryTreatmentInflow (const FEBasis& feBasis,
                        std::vector<bool>& dirichletNodes )
{
  static const int dim = FEBasis::GridView::dimension;

  // Interpolating the identity function wrt to a Lagrange basis
  // yields the positions of the Lagrange nodes
  Identity<dim> identity;
  BlockVector<FieldVector<double,dim> > lagrangeNodes;
  interpolate(feBasis, lagrangeNodes, identity);
  dirichletNodes.resize(lagrangeNodes.size());

  // Mark all Lagrange nodes on the bounding box as Dirichlet
  for (size_t i=0; i<lagrangeNodes.size(); i++)
  {
    bool isBoundary = false;
    for (int j=0; j<dim; j++)
      isBoundary = isBoundary || lagrangeNodes[i][j] < 1e-8; //|| lagrangeNodes[i][j] > bbox[j]-1e-8;

    if (isBoundary)

      dirichletNodes[i] = true;
  }
}


// A class implementing the analytical right hand side.  Here simply constant '1'
template <int dim>
class RightHandSide
    : public VirtualFunction<FieldVector<double,dim>, double >
{
public:
    void evaluate(const FieldVector<double,dim>& in, double& out) const {
        out = 1;
    }
};




int main(int argc, char** argv)
{
  try{
  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef YaspGrid<dim> GridType;
  FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {30, 30};
  GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  typedef Functions::Q2TraceNodalBasis<GridView> FEBasisTrace;              // u^
  FEBasisTrace feBasisTrace(gridView);

  typedef Functions::QKDiscontinuousNodalBasis<GridView, 1> FEBasisInterior;              // u
  FEBasisInterior feBasisInterior(gridView);

  typedef Functions::QKDiscontinuousNodalBasis<GridView, 4> FEBasisTest;              // v
  FEBasisTest feBasisTest(gridView);


  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////

  typedef BlockVector<FieldVector<double,1> > VectorType;
  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  VectorType rhs;
  MatrixType stiffnessMatrix;

  /////////////////////////////////////////////////////////
  //  Assemble the system
  /////////////////////////////////////////////////////////

  RightHandSide<dim> rightHandSide;
  assembleSystem<FEBasisTrace, FEBasisInterior, FEBasisTest>(feBasisTrace, feBasisInterior, feBasisTest, stiffnessMatrix, rhs, &rightHandSide);

  //printmatrix(std::cout , stiffnessMatrix, "stiffness", "--");

  /////////////////////////////////////////////////
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  VectorType x(feBasisTest.subIndexCount()+feBasisTrace.subIndexCount()+feBasisInterior.subIndexCount());
  x = 0;

  // Determine Dirichlet dofs for u^ (inflow boundary)
  std::vector<bool> dirichletNodesInflow;
  boundaryTreatmentInflow(feBasisTrace, dirichletNodesInflow);

/*  // Set Dirichlet values
  for (size_t i=0; i<feBasisTrace.subIndexCount(); i++)
  {
    if (dirichletNodesInflow[i])
    {
      // The zero is the value of the Dirichlet boundary condition
      rhs[feBasisTest.subIndexCount()+i] = 0;
    }
  }*/ //TODO unnoetig, da da eh schon 0 steht
  ////////////////////////////////////////////
  //   Modify Dirichlet rows
  ////////////////////////////////////////////

//printmatrix(std::cout , stiffnessMatrix, "stiffness1", "--");

  // loop over the matrix rows
  for (size_t i=0; i<feBasisTrace.subIndexCount(); i++)
  {
    if (dirichletNodesInflow[i])
    {
      auto cIt    = stiffnessMatrix[feBasisTest.subIndexCount()+i].begin();
      auto cEndIt = stiffnessMatrix[feBasisTest.subIndexCount()+i].end();
      // loop over nonzero matrix entries in current row
     for (; cIt!=cEndIt; ++cIt)
      {
        if (feBasisTest.subIndexCount()+i==cIt.index())
        {
          *cIt = 1;
        }
        else
        {
          *cIt = 0;
          stiffnessMatrix[cIt.index()][feBasisTest.subIndexCount()+i]=0;
        }
      }
    }

  }

/*// Determine Dirichlet dofs for v (inflow boundary)
  std::vector<bool> dirichletNodesInflowTest;
  boundaryTreatmentInflow(feBasisTest, dirichletNodesInflowTest);


  // Set Dirichlet values
  for (size_t i=0; i<feBasisTest.subIndexCount(); i++)
    if (dirichletNodesInflowTest[i])
      // The zero is the value of the Dirichlet boundary condition
      rhs[i] = 0;


  ////////////////////////////////////////////
  //   Modify Dirichlet rows
  ////////////////////////////////////////////

  // loop over the matrix rows
  for (size_t i=0; i<feBasisTest.subIndexCount(); i++)
  {
    if (dirichletNodesInflowTest[i])
    {
      auto cIt    = stiffnessMatrix[i].begin();
      auto cEndIt = stiffnessMatrix[i].end();
      // loop over nonzero matrix entries in current row
      for (; cIt!=cEndIt; ++cIt)
      {
        if (i==cIt.index())
        {
          *cIt = 1;
        }
        else
        {
          *cIt = 0;
          stiffnessMatrix[cIt.index()][i]=0;
        }
      }
    }
  }*/


/*  std::cout << "MatrixValues" << std::endl;
  for (unsigned int i=stiffnessMatrix.N()-4; i<stiffnessMatrix.N(); ++i)
  {
    for (unsigned int j=0; j<feBasisTest.subIndexCount(); ++j)
    {
      if(stiffnessMatrix.exists(i,j))
      {
        std::cout << std::scientific << std::setprecision(2) << stiffnessMatrix.entry(i,j) <<" ";
      }
      else
      {
        std::cout <<"         ";
      }
    }
    std::cout << std::endl;
  }*/
/*  std::cout << std::endl << "MatrixValues2" << std::endl;
    for (unsigned int i=0; i<stiffnessMatrix.N(); ++i)
    {
      for (unsigned int j=0; j<stiffnessMatrix.M(); ++j)
      {
        //if(std::abs(stiffnessMatrix.entry(i,j)-stiffnessMatrix.entry(j,i)) > 1e-4)
        if(stiffnessMatrix.exists(i,j))
        {
          if (stiffnessMatrix.entry(i,j) > 1e-4 || stiffnessMatrix.entry(i,j) < -1e-4)
          {
            std::cout <<"1";
          }
          else
          {
            std::cout <<" ";
          }
        }
        else
        {
          std::cout <<" ";
        }
      }
      std::cout << std::endl;
    }*/

  ////////////////////////////
  //   Compute solution
  ////////////////////////////

  std::cout <<"rhs size = "<< rhs.size()
            <<" matrix size = " << stiffnessMatrix.N() <<" x " << stiffnessMatrix.M()
            <<" solution size = "<< x.size() <<std::endl;



  //writeMatrixToMatlab(stiffnessMatrix, "TestMatrix1cell");

  UMFPack<MatrixType> umfPack(stiffnessMatrix, 2);
  InverseOperatorResult statistics;
  umfPack.apply(x, rhs, statistics);

#if 0
  // Technicality:  turn the matrix into a linear operator
  MatrixAdapter<MatrixType,VectorType,VectorType> op(stiffnessMatrix);


  // Sequential SSOR as preconditioner
  SeqSSOR<MatrixType,VectorType,VectorType> ssor(stiffnessMatrix, 3, 1);

  // Preconditioned BiCGSTABSolver
  BiCGSTABSolver<VectorType> BiCGSTAB(op,
                                      ssor, // preconditioner
                                      1e-4, // desired residual reduction factor
                                      50,   // maximum number of iterations
                                      2);   // verbosity of the solver)

  // Sequential incomplete LU decomposition as the preconditioner
/*  SeqILU0<MatrixType,VectorType,VectorType> ilu0(stiffnessMatrix,1.0);

  // Preconditioned conjugate-gradient solver
  CGSolver<VectorType> cg(op,
                          ilu0, // preconditioner
                          1e-4, // desired residual reduction factor
                          50,   // maximum number of iterations
                          2);   // verbosity of the solver*/

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  BiCGSTAB.apply(x, rhs, statistics);
#endif


  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////

//#if 0
  VectorType u(feBasisInterior.subIndexCount());
  u=0;
  for (unsigned int i=0; i<feBasisInterior.subIndexCount(); i++)
  {
    u[i] = x[feBasisTest.subIndexCount()+feBasisTrace.subIndexCount()+i];
  }

//  u=0;
  int idx = 2;
  //int idx =atoi(argv[1]);
//  u[idx]=1;

  //std::cout << u << std::endl;

  std::shared_ptr<VTKBasisGridFunction<FEBasisInterior,VectorType> > uFunction
    = std::make_shared<VTKBasisGridFunction<FEBasisInterior,VectorType> >(feBasisInterior, u, "solution");

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView,5);
  vtkWriter.addVertexData(uFunction);
  vtkWriter.write("solution_transport"+ std::to_string(idx));
//#endif
    return 0;
  }
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
