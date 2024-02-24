// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QKSUBSAMPLEDLOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_QKSUBSAMPLEDLOCALINTERPOLATION_HH

#include <dune/common/fvector.hh>
#include <dune/common/math.hh>
#include <dune/common/version.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>


namespace Dune
{
  /** \todo Please doc me! */
  template<int s, int k, int d, class LB>
  class QkSubsampledLocalInterpolation
  {

    // Return i as a d-digit number in the (s*k+1)-nary system
    static Dune::FieldVector<int,d> multiindex (int i)
    {
      Dune::FieldVector<int,d> alpha;
      for (int j=0; j<d; j++)
      {
        alpha[j] = i % (s*k+1);
        i = i/(s*k+1);
      }
      return alpha;
    }

  public:

    //! \brief Local interpolation of a function
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

      out.resize(power(s*k+1, d));

      for (int i = 0; i < power(s*k+1, d); i++)
      {
        // convert index i to multiindex
        Dune::FieldVector<int,d> alpha(multiindex(i));

        // Generate coordinate of the i-th Lagrange point
        for (int j=0; j<d; j++)
          x[j] = (1.0*alpha[j])/(s*k);

        out[i] = f(x);
      }
    }
  };

  /** \todo Please doc me! */
  template<int d, class LB>
  class QkSubsampledLocalInterpolation<1,0,d,LB>
  {
  public:
    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (
#if DUNE_VERSION_GTE(DUNE_LOCALFUNCTIONS,2,10)
        const F& f,
#else
        const F& ff,
#endif
        std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x(0);

#if DUNE_VERSION_LT(DUNE_LOCALFUNCTIONS,2,10)
      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);
#endif

      out.resize(1);
      out[0] = f(x);
    }
  };

}


#endif
