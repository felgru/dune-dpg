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

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#define DUNE_FINAL final

#include <dune/functions/functionspacebases/q2tracenodalbasis.hh>

#include <dune/dpg/functionhelper.hh>

using namespace Dune;







// Compute the stiffness matrix for a single element
template <class LocalView, class MatrixType>
void getLocalMatrix( const LocalView& localView, MatrixType& elementMatrix)
{
  // Get the grid element from the local FE basis view
  typedef typename LocalView::Element Element;
  const Element& element = localView.element();

  const int dim = Element::dimension;
  auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& localFiniteElement = localView.tree().finiteElement();

  // Set all matrix entries to zero
  elementMatrix.setSize(localFiniteElement.localBasis().size(),localFiniteElement.localBasis().size());
  elementMatrix = 0;      // fills the entire matrix with zeroes

  // Get a quadrature rule
  int order = 2*(dim-1);
  const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (size_t pt=0; pt < quad.size(); pt++) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The transposed inverse Jacobian of the map from the reference element to the element
    const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPos);

    // The gradients of the shape functions on the reference element
    std::vector<FieldMatrix<double,1,dim> > referenceGradients;
    localFiniteElement.localBasis().evaluateJacobian(quadPos, referenceGradients);

    // Compute the shape function gradients on the real element
    std::vector<FieldVector<double,dim> > gradients(referenceGradients.size());
    for (size_t i=0; i<gradients.size(); i++)
      jacobian.mv(referenceGradients[i][0], gradients[i]);

    // Compute the actual matrix entries
    for (size_t i=0; i<elementMatrix.N(); i++)
      for (size_t j=0; j<elementMatrix.M(); j++ )
        elementMatrix[i][j] += ( gradients[i] * gradients[j] ) * quad[pt].weight() * integrationElement;

  }

}


// Compute the source term for a single element
template <class LocalView>
void getVolumeTerm( const LocalView& localView,
                    BlockVector<FieldVector<double,1> >& localRhs,
                    const Dune::VirtualFunction<FieldVector<double,LocalView::Element::dimension>, double>* volumeTerm)
{
  // Get the grid element from the local FE basis view
  typedef typename LocalView::Element Element;
  const Element& element = localView.element();

  const int dim = Element::dimension;

  // Get set of shape functions for this element
  const auto& localFiniteElement = localView.tree().finiteElement();

  // Set all entries to zero
  localRhs.resize(localFiniteElement.localBasis().size());
  localRhs = 0;

  // A quadrature rule
  int order = dim;
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
    localFiniteElement.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

    // Actually compute the vector entries
    for (size_t i=0; i<localRhs.size(); i++)
      localRhs[i] += shapeFunctionValues[i] * functionValue * quad[pt].weight() * integrationElement;

  }

}

// Get the occupation pattern of the stiffness matrix
template <class InteriorTrialElement, class TestElement, class FEBasis>
void getOccupationPattern(const FEBasis& feBasis, MatrixIndexSet& nb)
{
  // Total number of grid vertices
  auto n = feBasis.subIndexCount();

  nb.resize(n, n);

  // A view on the FE basis on a single element
  typename FEBasis::LocalView localView(&feBasis);

  // Loop over all leaf elements
  auto it    = feBasis.gridView().template begin<0>();
  auto endIt = feBasis.gridView().template end<0>  ();

  for (; it!=endIt; ++it)
  {
    // Bind the local FE basis view to the current element
    localView.bind(*it);

        // There is a matrix entry a_ij if the i-th and j-th vertex are connected in the grid
        for (size_t i=0; i<localView.tree().localSize(); i++) {

            for (size_t j=0; j<localView.tree().localSize(); j++) {

                auto iIdx = localView.tree().globalIndex(i)[0];
                auto jIdx = localView.tree().globalIndex(j)[0];

                // Add a nonzero entry to the matrix
                nb.add(iIdx, jIdx);

            }

        }

    }

}


/** \brief Assemble the Laplace stiffness matrix on the given grid view */
template <class InteriorTrialElement, class TestElement, class FEBasis>
void assembleSystem(const FEBasis& feBasis,
                           BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                           BlockVector<FieldVector<double,1> >& rhs,
                           const VirtualFunction<Dune::FieldVector<double,FEBasis::GridView::dimension>,double>* volumeTerm)
{
  // Get the grid view from the finite element basis
  typedef typename FEBasis::GridView GridView;
  GridView gridView = feBasis.gridView();

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  MatrixIndexSet occupationPattern;
  getOccupationPattern<InteriorTrialElement, TestElement, FEBasis>(feBasis, occupationPattern);

  // ... and give it the occupation pattern we want.
  occupationPattern.exportIdx(matrix);

  // set rhs to correct length -- the total number of basis vectors in the basis
  rhs.resize(feBasis.subIndexCount());

  // Set all entries to zero
  matrix = 0;
  rhs = 0;

  // A view on the FE basis on a single element
  typename FEBasis::LocalView localView(&feBasis);

  // A loop over all elements of the grid
  auto it    = gridView.template begin<0>();
  auto endIt = gridView.template end<0>  ();

  for( ; it != endIt; ++it ) {

    // Bind the local FE basis view to the current element
    localView.bind(*it);

    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    getLocalMatrix(localView, elementMatrix);

    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i=0; i<elementMatrix.N(); i++) {

      // The global index of the i-th local degree of freedom of the element 'it'
      auto row = localView.tree().globalIndex(i)[0];

      for (size_t j=0; j<elementMatrix.M(); j++ ) {

        // The global index of the j-th local degree of freedom of the element 'it'
        auto col = localView.tree().globalIndex(j)[0];
        matrix[row][col] += elementMatrix[i][j];

      }

    }

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> > localRhs;
    getVolumeTerm(localView, localRhs, volumeTerm);

    for (size_t i=0; i<localRhs.size(); i++) {

      // The global index of the i-th vertex of the element 'it'
      auto row = localView.tree().globalIndex(i)[0];
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



// This method marks all vertices on the boundary of the grid.
// In our problem these are precisely the Dirichlet nodes.
// The result can be found in the 'dirichletNodes' variable.  There, a bit
// is set precisely when the corresponding vertex is on the grid boundary.
template <class FEBasis>
void boundaryTreatment (const FEBasis& feBasis,
                        const FieldVector<double,FEBasis::GridView::dimension>& bbox,
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
      isBoundary = isBoundary || lagrangeNodes[i][j] < 1e-8 || lagrangeNodes[i][j] > bbox[j]-1e-8;

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
  std::array<int,dim> elements = {2, 1};
  GridType grid(l,elements);

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid.leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  typedef Functions::Q2TraceNodalBasis<GridView> FEBasis;                    // u^
  FEBasis feBasis(gridView);

  typedef QkLocalFiniteElement<double,double,dim, 1> InteriorTrialElement;  // u
  typedef QkLocalFiniteElement<double,double,dim, 3> TestElement;           // v

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
  assembleSystem<InteriorTrialElement, TestElement, FEBasis>(feBasis, stiffnessMatrix, rhs, &rightHandSide);

  /////////////////////////////////////////////////
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  VectorType x(feBasis.subIndexCount());
  x = 0;

  // Determine Dirichlet dofs
  std::vector<bool> dirichletNodes;
  boundaryTreatment(feBasis, l, dirichletNodes);

  // Set Dirichlet values
  for (size_t i=0; i<rhs.size(); i++)
    if (dirichletNodes[i])
      // The zero is the value of the Dirichlet boundary condition
      rhs[i] = 0;

  ////////////////////////////////////////////
  //   Modify Dirichlet rows
  ////////////////////////////////////////////

  // loop over the matrix rows
  for (size_t i=0; i<stiffnessMatrix.N(); i++) {

    if (dirichletNodes[i]) {

      auto cIt    = stiffnessMatrix[i].begin();
      auto cEndIt = stiffnessMatrix[i].end();
      // loop over nonzero matrix entries in current row
      for (; cIt!=cEndIt; ++cIt)
        *cIt = (i==cIt.index()) ? 1 : 0;

    }

  }

  ////////////////////////////
  //   Compute solution
  ////////////////////////////

  // Technicality:  turn the matrix into a linear operator
  MatrixAdapter<MatrixType,VectorType,VectorType> op(stiffnessMatrix);

  // Sequential incomplete LU decomposition as the preconditioner
  SeqILU0<MatrixType,VectorType,VectorType> ilu0(stiffnessMatrix,1.0);

  // Preconditioned conjugate-gradient solver
  CGSolver<VectorType> cg(op,
                          ilu0, // preconditioner
                          1e-4, // desired residual reduction factor
                          50,   // maximum number of iterations
                          2);   // verbosity of the solver

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  cg.apply(x, rhs, statistics);

  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////

///  x=0;
//  int idx =atoi(argv[1]);
//  x[idx]=1;

  std::cout << x << std::endl;

  std::shared_ptr<VTKBasisGridFunction<FEBasis,VectorType> > xFunction
    = std::make_shared<VTKBasisGridFunction<FEBasis,VectorType> >(feBasis, x, "solution");

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView,5);
  vtkWriter.addVertexData(xFunction);
  vtkWriter.write("solution_transport");
    return 0;
  }
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
