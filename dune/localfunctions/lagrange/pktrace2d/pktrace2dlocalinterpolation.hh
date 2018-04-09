// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PKTRACE2DLOCALINTERPOLATION_HH
#define DUNE_PKTRACE2DLOCALINTERPOLATION_HH

#include <vector>

#include <dune/common/version.hh>

#if DUNE_VERSION_NEWER(DUNE_LOCALFUNCTIONS,2,7)
#include <dune/localfunctions/common/localinterpolation.hh>
#endif

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
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;

#if DUNE_VERSION_NEWER(DUNE_LOCALFUNCTIONS,2,7)
      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);
#else
      typename LB::Traits::RangeType y;
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
#if DUNE_VERSION_NEWER(DUNE_LOCALFUNCTIONS,2,7)
            out[n] = f(x);
#else
            ff.evaluate(x,y);
            out[n] = y;
#endif
            n++;
          }
        }
    }

  };
}

#endif
