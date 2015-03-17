#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <vector>

#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

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

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/dpg/assemble.hh>


using namespace Dune;



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

  // TODO: We are hacking our way around the fact that interpolation
  // of vector-value functions is not supported yet.
  BlockVector<FieldVector<double,dim> > lagrangeNodes;
  interpolate(feBasis, lagrangeNodes, [](FieldVector<double,dim> x){ return x; });

  dirichletNodes.resize(lagrangeNodes.size());

  // Mark all Lagrange nodes on the bounding box as Dirichlet
  for (size_t i=0; i<lagrangeNodes.size(); i++)
  {
    bool isBoundary = false;
    for (int j=0; j<dim; j++)
      isBoundary = isBoundary || lagrangeNodes[i][j] < 1e-8;

    if (isBoundary)

      dirichletNodes[i] = true;
  }
}




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

  typedef Functions::PQKTraceNodalBasis<GridView, 2> FEBasisTrace;            // u^
  FEBasisTrace feBasisTrace(gridView);

  typedef Functions::LagrangeDGBasis<GridView, 1> FEBasisInterior;              // u
  FEBasisInterior feBasisInterior(gridView);

  typedef Functions::LagrangeDGBasis<GridView, 4> FEBasisTest;              // v
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

  using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;  //TODO: WOZU?

  auto rightHandSide = [] (const Domain& x) { return 1;};
  //assembleSystem<FEBasisTrace, FEBasisInterior, FEBasisTest>(feBasisTrace, feBasisInterior, feBasisTest, stiffnessMatrix, rhs, rightHandSide);
  assembleSystem(feBasisTrace, feBasisInterior, feBasisTest, stiffnessMatrix, rhs, rightHandSide);

  /////////////////////////////////////////////////
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  VectorType x(feBasisTest.indexSet().size()+feBasisTrace.indexSet().size()
                +feBasisInterior.indexSet().size());
  x = 0;

  // Determine Dirichlet dofs for u^ (inflow boundary)
  std::vector<bool> dirichletNodesInflow;
  boundaryTreatmentInflow(feBasisTrace, dirichletNodesInflow);

  // TODO unnoetig, da da eh schon 0 steht
#if 0
  // Set Dirichlet values
  for (size_t i=0; i<feBasisTrace.indexSet().size(); i++)
  {
    if (dirichletNodesInflow[i])
    {
      // The zero is the value of the Dirichlet boundary condition
      rhs[feBasisTest.indexSet().size()+i] = 0;
    }
  }
#endif

  ////////////////////////////////////////////
  //   Modify Dirichlet rows
  ////////////////////////////////////////////

  // loop over the matrix rows
  for (size_t i=0; i<feBasisTrace.indexSet().size(); i++)
  {
    if (dirichletNodesInflow[i])
    {
      auto cIt    = stiffnessMatrix[feBasisTest.indexSet().size()+i].begin();
      auto cEndIt = stiffnessMatrix[feBasisTest.indexSet().size()+i].end();
      // loop over nonzero matrix entries in current row
     for (; cIt!=cEndIt; ++cIt)
      {
        if (feBasisTest.indexSet().size()+i==cIt.index())
        {
          *cIt = 1;
        }
        else
        {
          *cIt = 0;
          stiffnessMatrix[cIt.index()][feBasisTest.indexSet().size()+i]=0;
        }
      }
    }

  }

#if 0
  // Determine Dirichlet dofs for v (inflow boundary)
  std::vector<bool> dirichletNodesInflowTest;
  boundaryTreatmentInflow(feBasisTest, dirichletNodesInflowTest);


  // Set Dirichlet values
  for (size_t i=0; i<feBasisTest.indexSet().size(); i++)
    if (dirichletNodesInflowTest[i])
      // The zero is the value of the Dirichlet boundary condition
      rhs[i] = 0;


  ////////////////////////////////////////////
  //   Modify Dirichlet rows
  ////////////////////////////////////////////

  // loop over the matrix rows
  for (size_t i=0; i<feBasisTest.indexSet().size(); i++)
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
  }
#endif

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

  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////

  VectorType u(feBasisInterior.indexSet().size());
  u=0;
  for (unsigned int i=0; i<feBasisInterior.indexSet().size(); i++)
  {
    u[i] = x[feBasisTest.indexSet().size()+feBasisTrace.indexSet().size()+i];
  }

  Dune::Functions::DiscreteScalarGlobalBasisFunction<decltype(feBasisInterior),decltype(u)> uFunction(feBasisInterior,u);
  auto localUFunction = localFunction(uFunction);

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
  vtkWriter.addVertexData(localUFunction, VTK::FieldInfo("u", VTK::FieldInfo::Type::scalar, 1));
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
