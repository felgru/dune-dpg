// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PKFACE2DLOCALINTERPOLATION_HH
#define DUNE_PKFACE2DLOCALINTERPOLATION_HH

#include <vector>

#include <dune/common/version.hh>

#include <dune/localfunctions/common/localinterpolation.hh>

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
    void interpolate (
#if DUNE_VERSION_GTE(DUNE_LOCALFUNCTIONS,2,10)
        const F& f,
#else
        const F& ff,
#endif
        std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;

#if DUNE_VERSION_LT(DUNE_LOCALFUNCTIONS,2,10)
      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);
#endif

      out.resize(N);
      for (int l=0; l<=k; l++)
      {
        x[0]=1.0*l/kdiv;
        x[1]=0;
        out[l] = f(x);
        x[0]=0;
        x[1]=1.0*l/kdiv;
        out[l+(k+1)] = f(x);
        x[0]=1.0*l/kdiv;
        x[1]=1.0*(k-l)/kdiv;
        out[l+2*(k+1)] = f(x);
      }
    }

  };
}

#endif
