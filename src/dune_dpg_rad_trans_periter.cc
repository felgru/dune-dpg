#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <fstream>
#include <cstdlib> // for std::abort()

#include <vector>

#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/optimaltestbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/pqksubsampleddgbasis.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/dpg/radiative_transfer/periter.hh>
#include <dune/dpg/radiative_transfer/henyey_greenstein_scattering.hh>

#include <boost/math/constants/constants.hpp>


using namespace Dune;


// Value of the analytic solution "for the interior of the domain"
template <class Domain,class Direction>
double fInner(const Domain& x,
              const Direction& s)
{
  double c0 = 1.;
  double c1 = 1.;
  double value = (std::exp(c0*x[0])-1)*(std::exp(c1*x[1])-1);//v pure transport
  //double value = 1-(x[0]-0.5)*(x[0]-0.5)-(x[1]-0.5)*(x[1]-0.5); //v RT
  return value;
}
// Partial derivative of fInner with respect to x[0]
template <class Domain,class Direction>
double fInnerD0(const Domain& x,
                const Direction& s)
{
  double c0 = 1.;
  double c1 = 1.;
  double value = c0*std::exp(c0*x[0])*(std::exp(c1*x[1])-1);//v pure transport
  //double value = -2*(x[0]-0.5); //v RT
  return value;
}
// Partial derivative of fInner with respect to x[1]
template <class Domain,class Direction>
double fInnerD1(const Domain& x,
                const Direction& s)
{
  double c0 = 1.;
  double c1 = 1.;
  double value = (std::exp(c0*x[0])-1)*std::exp(c1*x[1])*c1;//v pure transport
  //double value = -2*(x[1]-0.5); //v RT
  return value;
}

// This function satifies the zero incoming flux bounday conditions
template <class Domain,class Direction>
double fBoundary(const Domain& x,
                 const Direction& s)
{
  // double value = 1.;//v pure transport
  double value = ( (s[0]>0)*x[0] + (s[0]==0)*1. + (s[0]<0)*(1-x[0]) )*
                 ( (s[1]>0)*x[1] + (s[1]==0)*1. + (s[1]<0)*(1-x[1]) );//v RT
  return value;
}
// Partial derivative of fBoundary with respect to x[0]
template <class Domain,class Direction>
double fBoundaryD0(const Domain& x,
                   const Direction& s)
{
  //double value = 0.;//v pure transport
  double value = ( (s[0]>0)*1 + (s[0]==0)*0. + (s[0]<0)*(-1.) )*
                 ( (s[1]>0)*x[1] + (s[1]==0)*1. + (s[1]<0)*(1-x[1]) );//v RT
  return value;
}
// Partial derivative of fBoundary with respect to x[1]
template <class Domain,class Direction>
double fBoundaryD1(const Domain& x,
                   const Direction& s)
{
  // double value = 0.;//v pure transport
  double value = ( (s[0]>0)*x[0] + (s[0]==0)*1. + (s[0]<0)*(1-x[0]) )*
                 ( (s[1]>0)*1 + (s[1]==0)*0. + (s[1]<0)*(-1.) ); //v RT
  return value;
}

//The analytic solution
template <class Domain, class Direction>
double uAnalytic(const Domain& x,
                 const Direction& s)
{
  return fInner(x,s)*fBoundary(x,s);
}

// Optical parameter: sigma
template <class Domain,class Direction>
double sigma(const Domain& x,
             const Direction& s)
{
  return 5.;
}

// Optical parameter: kernel
template <class Direction>
double kernel(const Direction& sIntegration,
              const Direction& s)
{
  return 1.7;
}

// The scattering kernel computed with an analytic solution u
// and a mid-point rule for the quadrature formula
template <class Domain, class Direction, class Function>
double collisionTerm(const Domain& x,
                     const Direction& s,
                     const Function& u,
                     const std::vector<Direction>& sVector)
{
  unsigned int numS = sVector.size();
  double sum = 0.;
  for(unsigned int i=0; i<numS; i++)
  {
    sum += kernel(sVector[i],s)*u(x,sVector[i]);
  }

  return sum/numS;
}

// The right hand-side
template <class Domain, class Direction, class Function>
double f(const Domain& x,
         const Direction& s,
         const Function& u,
         const std::vector<Direction>& sVector)
{
  double value = s[0]*( fInnerD0(x,s) * fBoundary(x,s) + fInner(x,s)*fBoundaryD0(x,s)) +
                 s[1]*( fInnerD1(x,s) * fBoundary(x,s) + fInner(x,s)*fBoundaryD1(x,s)) +
                 sigma(x,s)*u(x,s) - collisionTerm(x,s,u,sVector);
  return value;
}


// The values on the incoming boundary
template <class Domain,class Direction>
double gFunction(const Domain& x,
                 const Direction& s)
{
  return 1.5;
}


int main(int argc, char** argv)
{
  try{

  ///////////////////////////////////
  // Get arguments
  // argv[1]: number of discrete ordinates
  // argv[2]: number of fixed-point iterations
  // argv[3]: size of grid
  ///////////////////////////////////

  if(argc != 4) {
      std::cerr << "Usage: " << argv[0] << " <# of ordinates>"
                << " <# of iterations>"
                << " <size of grid>" << std::endl;
      std::abort();
  }

  unsigned int numS = atoi(argv[1]);
  int N = atoi(argv[2]);
  unsigned int sizeGrid = atoi(argv[3]);

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  const int dim = 2;
  typedef UGGrid<dim> GridType;

  FieldVector<double,dim> lower = {0,0};
  FieldVector<double,dim> upper = {1,1};
  array<unsigned int,dim> elements = {sizeGrid,sizeGrid};

  //shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, elements);

  shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createSimplexGrid(lower, upper, elements);

  //shared_ptr<GridType> grid = shared_ptr<GridType>(GmshReader<GridType>::read("irregular-square.msh"));

  using Domain = GridType::template Codim<0>::Geometry::GlobalCoordinate;
  using Direction = FieldVector<double, dim>;

  // Vector of directions: sVector
  std::vector<Direction> sVector(numS);
  for(unsigned int i = 0; i < numS; ++i)
  {
    using namespace boost::math::constants;
    sVector[i] = {cos(2*pi<double>()*i/numS),
                  sin(2*pi<double>()*i/numS)};
  }

  auto g = [&sVector] (const Domain& x, const Direction& s)
           { return 1.; };
  // TODO: Estimate ρ from the paper.
  const double rho = 1.;
  // TODO: Estimate the constant C_T.
  const double CT = 1;
  Periter<ScatteringKernelApproximation::HaarWavelet::MatrixCompression>()
      .solve(*grid, g, HenyeyGreensteinScattering<Direction>(0.9),
             numS, rho, CT, 1e-2, N);

  return 0;
  }
  catch (Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
