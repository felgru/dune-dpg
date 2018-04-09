// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PKFACE2DLOCALINTERPOLATION_HH
#define DUNE_PKFACE2DLOCALINTERPOLATION_HH

#include <vector>

#include <dune/common/version.hh>

#if DUNE_VERSION_NEWER(DUNE_LOCALFUNCTIONS,2,7)
#include <dune/localfunctions/common/localinterpolation.hh>
#endif

namespace Dune
{
  template<class LB>
  class PkFace2DLocalInterpolation
  {
    /** \brief The number of degrees of freedom */
    enum {N = LB::N};

    /** \brief Export the element order */
    enum {k = LB::O};

  private:
    enum {kdiv = (k == 0 ? 1 : k)};

  public:

    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;

#if DUNE_VERSION_NEWER(DUNE_LOCALFUNCTIONS,2,7)
      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);
#else
      typename LB::Traits::RangeType y;
#endif

      out.resize(N);
      for (int l=0; l<=k; l++)
      {
        x[0]=1.0*l/kdiv;
        x[1]=0;
#if DUNE_VERSION_NEWER(DUNE_LOCALFUNCTIONS,2,7)
        out[l] = f(x);
#else
        ff.evaluate(x,y);
        out[l] = y;
#endif
        x[0]=0;
        x[1]=1.0*l/kdiv;
#if DUNE_VERSION_NEWER(DUNE_LOCALFUNCTIONS,2,7)
        out[l+(k+1)] = f(x);
#else
        ff.evaluate(x,y);
        out[l+(k+1)] = y;
#endif
        x[0]=1.0*l/kdiv;
        x[1]=1.0*(k-l)/kdiv;
#if DUNE_VERSION_NEWER(DUNE_LOCALFUNCTIONS,2,7)
        out[l+2*(k+1)] = f(x);
#else
        ff.evaluate(x,y);
        out[l+2*(k+1)] = y;
#endif
      }
    }

  };
}

#endif
