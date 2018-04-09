// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <algorithm>
#include <limits>
#include <iostream>

#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/quadraturerules/subsampledquadraturerule.hh>

bool success = true;

/*
   This is a simple accuracy test on the reference element. It integrates
   x^p and y^p with the quadrature rule of order p, which should give
   an exact result.
 */

/*
   Exact (analytical) solution on different reference elements.
 */

template <class ctype, int dim>
ctype analyticalSolution (Dune::GeometryType t, int p, int direction )
{
  using Dune::GeometryType;
  ctype exact=0;

  if (t.isCube())
  {
    exact=1.0/(p+1);
    return exact;
  }

  if (t.isSimplex())
  {
    /* 1/(prod(k=1..dim,(p+k)) */
    exact = ctype( 1 );
    for( int k = 1; k <= dim; ++k )
      exact *= p+k;
    exact = ctype( 1 ) / exact;
    return exact;
  }

  if (t.isPrism())
  {
    const int pdim = (dim > 0 ? dim-1 : 0);
    if( direction < dim-1 )
    {
      GeometryType nt = Dune::GeometryTypes::simplex( dim-1 );
      if( dim > 0 )
        exact = analyticalSolution< ctype, pdim >( nt, p, direction );
      else
        exact = ctype( 1 );
    }
    else
      exact = ctype( 1 ) / ctype( Dune::Factorial< pdim >::factorial * (p+1));
    return exact;
  }

  if (t.isPyramid())
  {
    switch( direction )
    {
    case 0 :
    case 1 :
      exact=1.0/((p+3)*(p+1));
      break;
    case 2 :
      exact=2.0/((p+1)*(p+2)*(p+3));
      break;
    };
    return exact;
  }

  DUNE_THROW(Dune::NotImplemented, __func__ << " for " << t);
  return exact;
}

template<class QuadratureRule>
void checkQuadrature(const QuadratureRule &quad)
{
  using namespace Dune;
  typedef typename QuadratureRule::CoordType ctype;
  const unsigned int dim = QuadratureRule::d;
  const unsigned int p = quad.order();
  const Dune::GeometryType t = quad.type();
  FieldVector<ctype,dim> integral(0);
  for (const auto& qp : quad)
  {
    // pos of integration point
    const FieldVector< ctype, dim > &x = qp.position();
    const ctype weight = qp.weight();

    for (unsigned int d=0; d<dim; d++)
      integral[d] += weight*std::pow(x[d],double(p));
  }

  ctype maxRelativeError = 0;
  for(unsigned int d=0; d<dim; d++)
  {
    ctype exact = analyticalSolution<ctype,dim>(t,p,d);
    ctype relativeError = std::abs(integral[d]-exact) /
                          (std::abs(integral[d])+std::abs(exact));
    if (relativeError > maxRelativeError)
      maxRelativeError = relativeError;
  }
  ctype epsilon = std::pow(2.0,double(p))*p*std::numeric_limits<double>::epsilon();
  if (p==0)
    epsilon = 2.0*std::numeric_limits<double>::epsilon();
  if (maxRelativeError > epsilon) {
    std::cerr << "Error: Quadrature for " << t << " and order=" << p << " failed" << std::endl;
    for (unsigned int d=0; d<dim; d++)
    {
      ctype exact = analyticalSolution<ctype,dim>(t,p,d);
      ctype relativeError = std::abs(integral[d]-exact) /
                            (std::abs(integral[d])+std::abs(exact));
      std::cerr << "       relative error " << relativeError << " in direction " << d << " (exact = " << exact << " numerical = " << integral[d] << ")" << std::endl;
    }
    success = false;
  }
}

template<class QuadratureRule>
void checkWeights(const QuadratureRule &quad)
{
  typedef typename QuadratureRule::CoordType ctype;
  const unsigned int dim = QuadratureRule::d;
  const unsigned int p = quad.order();
  const Dune::GeometryType t = quad.type();
  double volume = 0;
  for (const auto& qp : quad)
  {
    volume += qp.weight();
  }
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
  if (std::abs(volume - Dune::referenceElement<ctype, dim>(t).volume())
      > quad.size()*std::numeric_limits<double>::epsilon())
#else
  if (std::abs(volume -
               Dune::ReferenceElements<ctype, dim>::general(t).volume())
      > quad.size()*std::numeric_limits<double>::epsilon())
#endif
  {
    std::cerr << "Error: Quadrature for " << t << " and order=" << p
              << " does not sum to volume of RefElem" << std::endl;
    std::cerr << "\tSums to " << volume << "( RefElem.volume() = "
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
              << Dune::referenceElement<ctype, dim>(t).volume()
              << ")" << "(difference "
              << volume - Dune::referenceElement<ctype, dim>(t).volume()
#else
              << Dune::ReferenceElements<ctype, dim>::general(t).volume()
              << ")" << "(difference " << volume -
                  Dune::ReferenceElements<ctype, dim>::general(t).volume()
#endif
              << ")" << std::endl;
    success = false;
  }
}

template<class ctype, int dim>
void check(Dune::GeometryType type,
           unsigned int maxOrder,
           Dune::QuadratureType::Enum qt = Dune::QuadratureType::GaussLegendre)
{
  typedef Dune::QuadratureRule<ctype, dim> Quad;

  for (unsigned int p=0; p<=maxOrder; ++p)
  {
    const Quad & quad = Dune::QuadratureRules<ctype,dim>::rule(type, p, qt);
    if (quad.type() != type || unsigned(quad.order()) < p) {
      std::cerr << "Error: Type mismatch! Requested Quadrature for " << type
                << " and order=" << p << "." << std::endl
                << "\tGot Quadrature for " << quad.type() << " and order="
                << quad.order() << std::endl;
      success = false;
      return;
    }
    checkWeights(quad);
    checkQuadrature(quad);
  }
  if (dim>0 && (dim>3 || type.isCube() || type.isSimplex()))
  {
    type = type.isCube() ? Dune::GeometryTypes::cube(dim-1) : Dune::GeometryTypes::simplex(dim-1);
    check<ctype,((dim==0) ? 0 : dim-1)>(type, maxOrder, qt);
  }
}

template<class ctype, int refinement, int dim>
void checkSubsampledRule(Dune::GeometryType type,
                         unsigned int maxOrder,
                         Dune::QuadratureType::Enum qt
                           = Dune::QuadratureType::GaussLegendre)
{
  typedef Dune::QuadratureRule<ctype, dim> BaseQuad;
  typedef Dune::SubsampledQuadratureRule<ctype, refinement, dim> Quad;

  for (unsigned int p=0; p<=maxOrder; ++p)
  {
    const BaseQuad& baseQuad = Dune::QuadratureRules<ctype,dim>::rule(type, p, qt);
    Quad quad = Quad(baseQuad);

    checkWeights(quad);
    checkQuadrature(quad);
  }
  if (dim>0 && (dim>3 || type.isCube() || type.isSimplex()))
  {
    type = type.isCube() ? Dune::GeometryTypes::cube(dim-1) : Dune::GeometryTypes::simplex(dim-1);
    check<ctype,((dim==0) ? 0 : dim-1)>(type, maxOrder, qt);
  }
}

int main (int argc, char** argv)
{
  unsigned int maxOrder = 45;
  if (argc > 1)
  {
    maxOrder = std::atoi(argv[1]);
    std::cout << "maxOrder = " << maxOrder << std::endl;
  }
  static const unsigned int refinement = 2;
  static const unsigned int dim = 2;

  checkSubsampledRule<double,refinement,dim>
      (Dune::GeometryTypes::simplex(dim), maxOrder);
  checkSubsampledRule<double,refinement,dim>
      (Dune::GeometryTypes::cube(dim), maxOrder);

  return success ? 0 : 1;
}
