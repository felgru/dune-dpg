#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#if ASTI_NORMALIZED_SPACES
#  define ASTI_SKELETAL_SCATTERING 1
#endif
#if ASTI_PEAKY_BV
#  define DUNE_DPG_USE_LEAST_SQUARES_INSTEAD_OF_CHOLESKY 1
#  define ASTI_SKELETAL_SCATTERING 1
#endif
#if ASTI_CHECKERBOARD
#  define ASTI_CHECKERBOARD_RHS 1
#  define ASTI_CHECKERBOARD_SIGMA 1
#endif
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>

#include <array>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/dpg/functions/gridviewfunctions.hh>
#ifndef ASTI_USE_UNIFORM_GRID
#  include <dune/dpg/radiative_transfer/asti.hh>
#else
#  include <dune/dpg/radiative_transfer/asti_uniform.hh>
#endif
#include <dune/dpg/radiative_transfer/henyey_greenstein_scattering.hh>


using namespace Dune;

void printHelp(const char* name) {
  std::cerr << "Usage: " << name
            << " [-psri] [-n <n>] [-o <dir>] [-k <k>] [-l <l>]"
            << " <target accuracy>"
            << " <gamma>"
            << " <# of iterations>"
            << " <size of grid>\n"
            << " -p: plot solutions\n"
            << " -s: plot scattering integral\n"
            << " -r: plot right hand side\n"
            << " -i: plot integrated solution\n"
            << " -n <n>: set maximal number of inner iterations to <n>\n"
            << " -o <dir>: set output directory to <dir>, default is "
               "\"./results/<current date and time>\"\n"
            << " -k <k>: set ratio between kernel and transport accuracy"
               " to <k>\n"
            << " -l <l>: set wavelet level of truth matrix for the scattering"
               " to <l>\n";
  std::exit(0);
}

void checkSizeGrid(unsigned int sizeGrid, unsigned int multipleOf) {
  if(sizeGrid == 0) {
    std::cerr << "Grid size has to be non-zero\n";
    std::exit(1);
  }
  if(sizeGrid % multipleOf != 0) {
    std::cerr << "Grid size has to be a multiple of " << multipleOf
              << " to resolve rhs and boundary data, but is "
              << sizeGrid << '\n';
    std::exit(1);
  }
}

int main(int argc, char** argv)
{
  // Set up MPI, if available
  FakeMPIHelper::instance(argc, argv);

  ///////////////////////////////////
  // Get arguments
  // argv[1]: target accuracy
  // argv[2]: gamma
  // argv[3]: maximal number of fixed-point iterations
  // argv[4]: size of grid
  ///////////////////////////////////

  ASTIPlotFlags plotFlags = ASTIPlotFlags::doNotPlot;
  std::string outdir = "";
  unsigned int maxNumberOfInnerIterations = 64;
  double accuracyRatio = 0.5;
  size_t scatteringTruthLevel = 6;

  {
    int opt;
    while ((opt = getopt(argc,argv,"n:o:k:l:psrih")) != EOF)
      switch(opt)
      {
        case 'p': plotFlags |= ASTIPlotFlags::plotOuterIterations; break;
        case 's': plotFlags |= ASTIPlotFlags::plotScattering; break;
        case 'r': plotFlags |= ASTIPlotFlags::plotRhs; break;
        case 'i': plotFlags |= ASTIPlotFlags::plotIntegratedSolution; break;
        case 'n': maxNumberOfInnerIterations = std::atoi(optarg); break;
        case 'o': outdir = optarg; break;
        case 'k': accuracyRatio = std::atof(optarg); break;
        case 'l': scatteringTruthLevel = std::atoi(optarg); break;
        default:
        case '?':
        case 'h':
          printHelp(argv[0]);
      }
    if(optind != argc-4) {
      printHelp(argv[0]);
    }
  }

  const unsigned int wltOrder = 2;
  const double targetAccuracy = std::atof(argv[optind]);
  const double gamma = std::atof(argv[optind+1]);
  const unsigned int N = std::atoi(argv[optind+2]);
  const unsigned int sizeGrid = std::atoi(argv[optind+3]);
#if ASTI_PEAKY_BV
  checkSizeGrid(sizeGrid, 8);
#endif
#if ASTI_CHECKERBOARD
  checkSizeGrid(sizeGrid, 7);
#endif

  if (outdir == "") {
    const auto now = std::chrono::system_clock::now();
    const std::time_t cnow = std::chrono::system_clock::to_time_t(now);
    std::stringstream folderstream;
    folderstream << "results/"
                 << std::put_time(std::localtime(&cnow), "%F-time%H%M%S");
    outdir = folderstream.str();
  }
  std::system(("mkdir -p "+outdir).data());
  std::cout << "output path: " << outdir << "\n\n";

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  constexpr int dim = 2;
  using Grid = UGGrid<dim>;
  using GridView = typename Grid::LeafGridView;

  using Domain = Grid::template Codim<0>::Geometry::GlobalCoordinate;
  using Direction = FieldVector<double, dim>;

  const FieldVector<double,dim> lower = {0,0};
  const FieldVector<double,dim> upper = {1,1};
  const std::array<unsigned int,dim> elements = {sizeGrid,sizeGrid};

  //std::unique_ptr<Grid> grid = StructuredGridFactory<Grid>::createCubeGrid(lower, upper, elements);

  std::unique_ptr<Grid> grid = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, elements);

  //std::unique_ptr<Grid> grid = GmshReader<Grid>::read("irregular-square.msh");

#ifndef ASTI_USE_UNIFORM_GRID
  // UG by default uses red-green refinements which would create and remove
  // auxiliary cells. This doesn't play well with the SubGrids we use, so
  // disable it here.
  grid->setClosureType(Grid::NONE);
#endif

#if ASTI_CHECKERBOARD_RHS
  const auto f
    = [](const GridView gridView)
      {
        auto rhs = [](const Domain& x)
        {
          constexpr int n=7;
          const int i = static_cast<int>(n*x[0]);
          const int j = static_cast<int>(n*x[1]);
          constexpr double v1 = 0.;
          constexpr double v2 = 1.;
          if(i==3 and j==3) {
            return v2;
          } else {
            return v1;
          }
        };
        return Functions::makePiecewiseConstantGridViewFunction(
              std::move(rhs), gridView);
      };
#else
  const auto f
    = [](const GridView)
      { return [](const Domain& x){ return 1.; }; };
#endif

  using RHSApproximation = FeRHS;

#if ASTI_CHECKERBOARD_SIGMA
  const auto sigma = [](const auto gridView)
      {
        auto sigma_ = [](const Domain& x)
            {
              constexpr int n=7;
              const int i = static_cast<int>(n*x[0]);
              const int j = static_cast<int>(n*x[1]);
              constexpr double v1 = 2.;
              constexpr double v2 = 10.;
              if(i<=0 or i>=6 or j<=0 or j>=6 or
                  (i+j) % 2 != 0 or (i==3 and (j==3 or j==5))) {
                return v1;
              } else {
                return v2;
              }
            };
        return Functions::makePiecewiseConstantGridViewFunction(
              std::move(sigma_), gridView);
      };
  const double sigmaMin = 2.;
  const double sigmaMax = 10.;
#else
  constexpr double sigma_ = 5.;
  const auto sigma = [sigma_](const auto gridView)
      {
        return Functions::makeConstantGridViewFunction(sigma_, gridView);
      };
  const double sigmaMin = sigma_;
  const double sigmaMax = sigma_;
#endif

#if ASTI_PEAKY_BV
  const auto g = [](const Domain& x)
      {
        if(x[0] < .1) {
          const double xProj = x[1];
          if(xProj>=.5+.125 || xProj<=.5-.125) {
            return 0.;
          } else if(xProj<=.5) {
            return 8*(xProj-(.5-.125));
          } else {
            return 1-8*(xProj-.5);
          }
        } else {
          return 0.;
        }
      };
  const auto homogeneous_inflow_boundary =
    [](const Direction& s) { return !(s[0] > 0.); };
#else
  const auto g = [](const Domain& x) { return 0.; };
  const auto homogeneous_inflow_boundary =
    [](const Direction& s) { return true; };
#endif

  const double domainDiameter = std::sqrt(2.);
  // Formula from Lemma 2.8 paper [DGM]
  const double CT
    = std::min(domainDiameter, std::sqrt(domainDiameter/(2*sigmaMin)));
  // Formula from Lemma 2.13 paper [DGM]
  const double rho
    = std::min(sigmaMax/(sigmaMin*sigmaMin),
        std::min((sigmaMax-sigmaMin+1.)/sigmaMin,
          std::sqrt(domainDiameter/(2*sigmaMin))));
  assert(rho < 1.);
  // using Proposition 2.11 from our paper [DGM]
  const double cB = sigmaMin - 1.;

  ///////////////////////////////////
  // Parameters for adaptivity
  ///////////////////////////////////

  // TODO: estimate norm of rhs f in V'
  // Remark: Here, V=H_{0,+}(D\times S)
#if ASTI_CHECKERBOARD_RHS
  const double fnorm = 1./7.;
#else
  const double fnorm = 1;
#endif
  const double err0 = fnorm / cB;

  ASTIApproximationParameters approximationParameters(accuracyRatio,
                                                         rho, CT, err0,
                                                         scatteringTruthLevel,
                                                         RHSApproximation{});

  ASTI<ScatteringKernelApproximation::AlpertWavelet::SVD<wltOrder>,
          RHSApproximation>()
      .solve(*grid, f, g, homogeneous_inflow_boundary, sigma,
             HenyeyGreensteinScattering(gamma),
             approximationParameters, targetAccuracy, N,
             maxNumberOfInnerIterations, outdir, plotFlags);

  return 0;
}
