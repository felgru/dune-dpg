// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PKFACE2DLOCALINTERPOLATION_HH
#define DUNE_PKFACE2DLOCALINTERPOLATION_HH

#include <vector>

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
    static const int kdiv = (k == 0 ? 1 : k);

  public:

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;
      typedef typename LB::Traits::DomainFieldType D;
      out.resize(N);
      int n=0;
      for (int l=0; l<=k; l++)
      {
        x[0]=1.0*l/kdiv;
        x[1]=0;
        f.evaluate(x,y);
        out[l] = y;
        x[0]=0;
        x[1]=1.0*l/kdiv;
        f.evaluate(x,y);
        out[l+(k+1)] = y;
        x[0]=1.0*l/kdiv;
        x[1]=1.0*(k-l)/kdiv;
        f.evaluate(x,y);
        out[l+2*(k+1)] = y;
      }
    }

  };
}

#endif
