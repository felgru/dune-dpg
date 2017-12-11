#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <algorithm>
#include <chrono>
#include <cstdlib> // for std::exit() and std::system()
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>

#include <array>
#include <vector>

#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#ifndef PERITER_USE_UNIFORM_GRID
#  include <dune/dpg/radiative_transfer/periter.hh>
#else
#  include <dune/dpg/radiative_transfer/periter_uniform.hh>
#endif
#include <dune/dpg/radiative_transfer/henyey_greenstein_scattering.hh>


using namespace Dune;

void printHelp(const char* name) {
  std::cerr << "Usage: " << name
            << " [-p] [-s] [-n <n>] [-o <dir>]"
            << " <target accuracy>"
            << " <gamma>"
            << " <# of iterations>"
            << " <size of grid>\n"
            << " -p: plot solutions\n"
            << " -s: plot scattering integral\n"
            << " -n <n>: set maximal number of inner iterations to <n>\n"
            << " -o <dir>: set output directory to <dir>, default is "
               "\"../results/\"\n";
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
  try{

  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  ///////////////////////////////////
  // Get arguments
  // argv[1]: target accuracy
  // argv[2]: gamma
  // argv[3]: maximal number of fixed-point iterations
  // argv[4]: size of grid
  ///////////////////////////////////

  PlotSolutions plotSolutions = PlotSolutions::doNotPlot;
  std::string basedir = "../results/";
  unsigned int maxNumberOfInnerIterations = 64;

  {
    int opt;
    while ((opt = getopt(argc,argv,"n:o:psh")) != EOF)
      switch(opt)
      {
        case 'p': plotSolutions |= PlotSolutions::plotOuterIterations; break;
        case 's': plotSolutions |= PlotSolutions::plotScattering; break;
        case 'n': maxNumberOfInnerIterations = atoi(optarg); break;
        case 'o': basedir = optarg; break;
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
  const double targetAccuracy = atof(argv[optind]);
  const double gamma = atof(argv[optind+1]);
  const unsigned int N = atoi(argv[optind+2]);
  const unsigned int sizeGrid = atoi(argv[optind+3]);
#if PERITER_PEAKY_BV
  checkSizeGrid(sizeGrid, 8);
#elif PERITER_CHECKERBOARD
  checkSizeGrid(sizeGrid, 7);
#endif

  std::string foldername;
  {
    const std::chrono::system_clock::time_point now
        = std::chrono::system_clock::now();
    const std::time_t cnow = std::chrono::system_clock::to_time_t(now);
    std::stringstream folderstream;
    folderstream << basedir
                 << std::put_time(std::localtime(&cnow), "%F-time%H%M%S");
    foldername = folderstream.str();
  }
  std::system(("mkdir -p "+foldername).data());
  std::cout << "output path: " << foldername << "\n\n";

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef UGGrid<dim> GridType;

  using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;
  using Direction = FieldVector<double, dim>;

  const FieldVector<double,dim> lower = {0,0};
  const FieldVector<double,dim> upper = {1,1};
  const std::array<unsigned int,dim> elements = {sizeGrid,sizeGrid};

  //std::shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  std::shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);

  //std::shared_ptr<GridType> grid = std::shared_ptr<GridType>(GmshReader<GridType>::read("irregular-square.msh"));

#ifndef PERITER_USE_UNIFORM_GRID
  // UG by default uses red-green refinements which would create and remove
  // auxiliary cells. This doesn't play well with the SubGrids we use, so
  // disable it here.
  grid->setClosureType(GridType::NONE);
#endif

#if PERITER_PEAKY_BV
  const auto f
    = [](const Domain& x, const Direction& s)
      { return 1.; };
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
  const double sigma = 5.;
  const double sigmaMin = sigma;
  const double sigmaMax = sigma;
#elif PERITER_CHECKERBOARD
  const auto f
    = [](const Domain& x, const Direction& s)
      {
        const int n=7;
        const int i = static_cast<int>(n*x[0]);
        const int j = static_cast<int>(n*x[1]);
        const double v1 = 0.;
        const double v2 = 1.;
        if(i==3 and j==3) {
          return v2;
        } else {
          return v1;
        }
      };
  const auto g = [](const Domain& x) { return 0.; };
  const auto homogeneous_inflow_boundary =
    [](const Direction& s) { return true; };
  const auto sigma = [](const Domain& x)
      {
        const int n=7;
        const int i = static_cast<int>(n*x[0]);
        const int j = static_cast<int>(n*x[1]);
        const double v1 = 2.;
        const double v2 = 10.;
        if(i<=0 or i>=6 or j<=0 or j>=6 or
            (i+j) % 2 == 0 or (i==3 and j==5)) {
          return v1;
        } else {
          return v2;
        }
      };
  const double sigmaMin = 2.;
  const double sigmaMax = 10.;
#else
#  error "Not specified which problem to solve."
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

#if PERITER_PEAKY_BV
#  define DUNE_DPG_USE_LEAST_SQUARES_INSTEAD_OF_CHOLESKY 1
#endif

  Periter<ScatteringKernelApproximation::AlpertWavelet::SVD<wltOrder>,
          FeRHS>()
      .solve(*grid, f, g, homogeneous_inflow_boundary, sigma,
             HenyeyGreensteinScattering(gamma),
             rho, CT, cB, targetAccuracy, N, maxNumberOfInnerIterations,
             foldername, plotSolutions);

  return 0;
  }
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
