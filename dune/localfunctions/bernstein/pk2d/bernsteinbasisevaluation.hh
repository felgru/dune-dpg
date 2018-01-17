// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BERNSTEIN_PK2D_BASIS_EVALUTATION_HH
#define DUNE_BERNSTEIN_PK2D_BASIS_EVALUTATION_HH

#include <array>
#include <functional>
#include <numeric>

#include <dune/common/fvector.hh>

namespace Dune {
namespace detail {

  template<class D, class R, unsigned int k>
  struct Pk2DBernsteinBasis
  {
    using BarycentricCoordinate = Dune::FieldVector<D,3>;

    template<class DomainType, class OutputIterator>
    static void fillVectorOfEvaluations
        (const DomainType& x, OutputIterator&& out)
    {
      const BarycentricCoordinate xBarycentric{x[0], x[1], 1-x[0]-x[1]};
      BarycentricCoordinate xPowers{1., 1., 1.};
      const auto factorial = factorials();
      for (unsigned int i0=0; i0<=k; i0++, xPowers[0] *= xBarycentric[0])
      {
        xPowers[1] = 1.;
        for (unsigned int i1=0; i1<=k-i0; i1++, xPowers[1] *= xBarycentric[1])
        {
          const unsigned int i2 = k-i0-i1;
          xPowers[2] = 1.;
          for(unsigned int l=0; l<i2; l++)
            xPowers[2] *= xBarycentric[2];
          *out = std::accumulate(xPowers.begin(), xPowers.end(),
                    static_cast<R>(factorial[k])
                      / (factorial[i0] * factorial[i1] * factorial[i2]),
                    std::multiplies<R>());
          ++out;
        }
      }
    }

    static std::array<unsigned int, k+1> factorials()
    {
      std::array<unsigned int, k+1> factorials;
      factorials[0] = 1;
      unsigned int last = 1;
      for(unsigned int l=1; l<=k; l++) {
        last *= l;
        factorials[l] = last;
      }
      return factorials;
    }
  };

  // specialization for k==0, not clear whether that is needed
  template<class D, class R>
  struct Pk2DBernsteinBasis<D, R, 0>
  {
    template<class DomainType, class OutputIterator>
    static void fillVectorOfEvaluations
        (const DomainType& x, OutputIterator&& out)
    {
      *out = 1;
    }
  };

  // TODO: we might need a specialization for k==-1, but k is unsigned
  template<class D, class R>
  struct Pk2DBernsteinBasis<D, R, -1>
  {
    template<class DomainType, class OutputIterator>
    static void fillVectorOfEvaluations
        (const DomainType& x, OutputIterator&& out)
    {
      *out = 0;
    }
  };

}} // end namespace Dune::detail
#endif
