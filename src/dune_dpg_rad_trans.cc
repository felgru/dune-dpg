#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <cstdlib> // for std::abort()

#include <vector>

#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/common/intersection.hh> //TODO necessary?
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/gmshreader.hh>



#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/istl/umfpack.hh>


#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/optimaltestbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/dpg/system_assembler.hh>
#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/radiative_transfer/scattering.hh>

#include <boost/math/constants/constants.hpp>


using namespace Dune;



int main(int argc, char** argv)
{
  try{
  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef UGGrid<dim> GridType;

  FieldVector<double,dim> lower = {0,0};
  FieldVector<double,dim> upper = {1,1};
  array<unsigned int,dim> elements = {10,10};

  //shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);

  //shared_ptr<GridType> grid = shared_ptr<GridType>(GmshReader<GridType>::read("irregular-square.msh"));

  typedef GridType::LeafGridView GridView;
  GridView gridView = grid->leafGridView();

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  typedef Functions::LagrangeDGBasis<GridView, 1> FEBasisInterior; // u
  FEBasisInterior feBasisInterior(gridView);

  typedef Functions::PQKTraceNodalBasis<GridView, 2> FEBasisTrace; // u^
  FEBasisTrace feBasisTrace(gridView);

  auto solutionSpaces = std::make_tuple(FEBasisInterior(gridView), FEBasisTrace(gridView));

  typedef Functions::LagrangeDGBasis<GridView, 4> FEBasisTest;     // v enriched
  auto testSpaces = std::make_tuple(FEBasisTest(gridView));

  typedef decltype(testSpaces) TestSpaces;
  typedef decltype(solutionSpaces) SolutionSpaces;

  typedef decltype(make_BilinearForm(testSpaces, solutionSpaces,
            make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(0.),
              make_IntegralTerm<0,0,IntegrationType::gradValue,
                                    DomainOfIntegration::interior>(-1.,
                                       FieldVector<double, dim>{1.,1.}),
              make_IntegralTerm<0,1,IntegrationType::normalVector,
                                    DomainOfIntegration::face>(1.,
                                       FieldVector<double, dim>{1.,1.}))))
          BilinearForm;
  typedef decltype(make_InnerProduct(testSpaces,
            make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(1.),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>(1.,
                                       FieldVector<double, dim>{1.,1.}))))
          InnerProduct;

  typedef Functions::TestspaceCoefficientMatrix<BilinearForm, InnerProduct>
      TestspaceCoefficientMatrix;
  typedef Functions::OptimalTestBasis<TestspaceCoefficientMatrix>
      FEBasisOptimalTest;              // v

  typedef decltype(make_SystemAssembler(
#if 1
              std::declval<std::tuple<FEBasisOptimalTest>>(), solutionSpaces,
#else
              testSpaces, solutionSpaces,
#endif
              std::declval<BilinearForm>(),
              std::declval<InnerProduct>(), DPGFormulation()))
          SystemAssembler_t;

  if(argc != 3) {
      std::cerr << "Usage: " << argv[0] << " <# of ordinates>"
                << " <# of iterations>" << std::endl;
      std::abort();
  }
  int numS = atoi(argv[1]);
  int N = atoi(argv[2]);

  std::vector<SystemAssembler_t> systemAssemblers;
  systemAssemblers.reserve(numS);

  std::vector<ScatteringAssembler<std::tuple<FEBasisOptimalTest>,
                                  SolutionSpaces,
                                  DPGFormulation>
             > scatteringAssemblers;
  scatteringAssemblers.reserve(numS);

  /* create an FEBasisOptimalTest for each direction */
  std::vector<std::tuple<FEBasisOptimalTest> > optimalTestSpaces;
  optimalTestSpaces.reserve(numS);
  /* All the following objects have to be created outside of the
   * following for loop, as the optimalTestSpace holds references
   * to them which will otherwise go out of scope. */
  std::vector<BilinearForm> bilinearForms;
  bilinearForms.reserve(numS);
  std::vector<InnerProduct> innerProducts;
  innerProducts.reserve(numS);
  std::vector<TestspaceCoefficientMatrix> coefficientMatrices;
  coefficientMatrices.reserve(numS);

  for(int i = 0; i < numS; ++i)
  {
    using namespace boost::math::constants;
    FieldVector<double, dim> s = {cos(2*pi<double>()*i/numS),
                                  sin(2*pi<double>()*i/numS)};
    bilinearForms.emplace_back(
      make_BilinearForm(testSpaces, solutionSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(2.),
              make_IntegralTerm<0,0,IntegrationType::gradValue,
                                    DomainOfIntegration::interior>(-1., s),
              make_IntegralTerm<0,1,IntegrationType::normalVector,
                                    DomainOfIntegration::face>(1., s))));
    innerProducts.emplace_back(
      make_InnerProduct(testSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::valueValue,
                                    DomainOfIntegration::interior>(1.),
              make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>(1., s))));

    coefficientMatrices.emplace_back(bilinearForms[i], innerProducts[i]);

    optimalTestSpaces.emplace_back(
            make_tuple(FEBasisOptimalTest(coefficientMatrices[i])));

    systemAssemblers.emplace_back(
        make_SystemAssembler(optimalTestSpaces[i], solutionSpaces,
                             bilinearForms[i], innerProducts[i],
                             DPGFormulation()));
    scatteringAssemblers.emplace_back(
        make_ScatteringAssembler(optimalTestSpaces[i],
                                 solutionSpaces,
                                 DPGFormulation()));
  }

  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////


  typedef BlockVector<FieldVector<double,1> > VectorType;
  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  std::vector<VectorType> rhs(numS);
  std::vector<MatrixType> stiffnessMatrix(numS);

  // Determine Dirichlet dofs for u^ (inflow boundary)
  std::vector<std::vector<bool>> dirichletNodesInflow(numS);
  for(int i = 0; i < numS; ++i)
  {
    using namespace boost::math::constants;
    FieldVector<double, dim> s = {cos(2*pi<double>()*i/numS),
                                  sin(2*pi<double>()*i/numS)};
    BoundaryTools boundaryTools = BoundaryTools();
    boundaryTools.boundaryTreatmentInflow(std::get<1>(solutionSpaces),
                                          dirichletNodesInflow[i],
                                          s);
  }

  /////////////////////////////////////////////////
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  std::vector<VectorType> x;
  x.reserve(numS);
  for(int i = 0; i < numS; ++i)
  {
    x.emplace_back(feBasisTrace.indexSet().size()
                   +feBasisInterior.indexSet().size());
    x[i] = 0;
  }

  for(int n = 0; n < N; ++n)
  {

    /////////////////////////////////////////////////////////
    //  Assemble the systems
    /////////////////////////////////////////////////////////
    using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;
    auto f = std::make_tuple([] (const Domain& x) { return 1.;});
    for(int i = 0; i < numS; ++i)
    {
      systemAssemblers[i].assembleSystem(stiffnessMatrix[i], rhs[i], f);
      VectorType scattering;
      scatteringAssemblers[i].assembleScattering<0>(scattering, i, x);
      rhs[i] += scattering;
      systemAssemblers[i].applyDirichletBoundarySolution<1>
          (stiffnessMatrix[i],
           rhs[i],
           dirichletNodesInflow[i],
           0.);
    }

    // std::ofstream of("stiffnessNew.dat");
    // printmatrix(of, stiffnessMatrix[0], "stiffnessNew", "--");

    ////////////////////////////
    //   Compute solution
    ////////////////////////////

    std::cout <<"rhs size = "<< rhs[0].size()
              <<" matrix size = " << stiffnessMatrix[0].N()
                         << " x " << stiffnessMatrix[0].M()
              <<" solution size = "<< x[0].size() <<std::endl;


    for(int i = 0; i < numS; ++i)
    {
      UMFPack<MatrixType> umfPack(stiffnessMatrix[i], 2);
      InverseOperatorResult statistics;
      umfPack.apply(x[i], rhs[i], statistics);
    }


    //////////////////////////////////////////////////////////////////////////
    //  Make a discrete function from the FE basis and the coefficient vector
    //////////////////////////////////////////////////////////////////////////

    for(int i = 0; i < numS; ++i)
    {
      VectorType u(feBasisInterior.indexSet().size());
      u=0;
      for (unsigned int j=0; j<feBasisInterior.indexSet().size(); j++)
      {
        u[j] = x[i][j];
      }

      Dune::Functions::DiscreteScalarGlobalBasisFunction
          <decltype(feBasisInterior),decltype(u)>
          uFunction(feBasisInterior,u);
      auto localUFunction = localFunction(uFunction);

      ////////////////////////////////////////////////////////////////////////
      //  Write result to VTK file
      //  We need to subsample, because VTK cannot natively display
      //  real second-order functions
      ////////////////////////////////////////////////////////////////////////
      SubsamplingVTKWriter<GridView> vtkWriter(gridView,2);
      vtkWriter.addVertexData(localUFunction,
                      VTK::FieldInfo("u", VTK::FieldInfo::Type::scalar, 1));
      std::string name = std::string("solution_rad_trans_n")
                       + std::to_string(n)
                       + std::string("_s")
                       + std::to_string(i);
      vtkWriter.write(name);

    }
  }

  return 0;
  }
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
