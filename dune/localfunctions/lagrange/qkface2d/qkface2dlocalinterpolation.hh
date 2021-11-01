// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QKFACE2DLOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_QKFACE2DLOCALINTERPOLATION_HH

#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>


namespace Dune
{
  /** \todo Please doc me! */
  template<int k, class LB>
  class QkFace2DLocalInterpolation
  {

  public:

    //! \brief Local interpolation of a function -> works only for functions with support only on the boundary
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;

      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

      out.resize(4*(k+1));
      for (unsigned int l=0; l<=k; l++)
      {
        x[0]=0;
        x[1]=1.0*l/k;
        out[l] = f(x);
        x[0]=1;
        x[1]=1.0*l/k;
        out[l+(k+1)] = f(x);
        x[0]=1.0*l/k;
        x[1]=0;
        out[l+2*(k+1)] = f(x);
        x[0]=1.0*l/k;
        x[1]=1;
        out[l+3*(k+1)] = f(x);
      }
    }
  };

  /** \todo Please doc me! */
  template<class LB>
  class QkFace2DLocalInterpolation<0,LB>
  {
  public:
    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x(0);

      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

      out.resize(1);
      out[0] = f(x);
    }
  };
}


#endif
