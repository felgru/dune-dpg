// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PKTRACE2DLOCALINTERPOLATION_HH
#define DUNE_PKTRACE2DLOCALINTERPOLATION_HH

#include <vector>

#include <dune/common/version.hh>

#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{
  template<class LB>
  class PkTrace2DLocalInterpolation
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

      typedef typename LB::Traits::DomainFieldType D;
      out.resize(N);
      int n=0;
      for (int j=0; j<=k; j++)
        for (int i=0; i<=k-j; i++)
        {
          if (i==0 or i+j==k or j==0)
          {
            x[0] = static_cast<D>(i)/static_cast<D>(kdiv);
            x[1] = static_cast<D>(j)/static_cast<D>(kdiv);
            out[n] = f(x);
            n++;
          }
        }
    }

  };
}

#endif
