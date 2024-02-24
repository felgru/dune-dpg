// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QKTRACELOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_QKTRACELOCALINTERPOLATION_HH

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
  template<int k, int d, class LB>
  class QkTraceLocalInterpolation
  {

    // Return i as a d-digit number in the (k+1)-nary system
    static Dune::FieldVector<int,d> multiindex (int i)
    {
      Dune::FieldVector<int,d> alpha;
      for (int j=0; j<d; j++)
      {
        alpha[j] = i % (k+1);
        i = i/(k+1);
      }
      return alpha;
    }

  public:

    //! \brief Local interpolation of a function -> works only for functions with support only on the boundary
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

      out.resize(power(k+1, d) - power(k-1, d));
      unsigned int i = 0;
      for (int l = 0; l < power(k+1, d); l++)
      {
        // convert index i to multiindex
        const Dune::FieldVector<int,d> alpha(multiindex(l));
       // check if the corresponding lagrange polynomial has support on any face
        bool is_on_face=false;
        for (unsigned int b=0; b<d && !is_on_face; b++)
        {
          is_on_face=(alpha[b]==0 || alpha[b]==k);
        }
        if (is_on_face)
        {
          // Generate coordinate of the i-th Lagrange point
          for (int j=0; j<d; j++)
          {
            x[j] = (1.0*alpha[j])/k;
          }
          out[i] = f(x);
          i++;
        }
      }
    }
  };

  /** \todo Please doc me! */
  template<int d, class LB>
  class QkTraceLocalInterpolation<0,d,LB>
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
