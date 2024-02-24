// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BERNSTEIN_PK2DLOCALINTERPOLATION_HH
#define DUNE_BERNSTEIN_PK2DLOCALINTERPOLATION_HH

#include <algorithm>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <dune/localfunctions/common/localinterpolation.hh>

#include "bernsteinbasisevaluation.hh"

namespace Dune
{
  template<class LB>
  class BernsteinPk2DLocalInterpolation
  {
    /** \brief The number of degrees of freedom */
    enum {N = LB::N};

    /** \brief Export the element order */
    enum {k = LB::O};

  private:
    static const int kdiv = (k == 0 ? 1 : k);

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
      typedef typename LB::Traits::DomainFieldType D;
      typedef typename LB::Traits::RangeFieldType R;

      typename LB::Traits::DomainType x;
#if DUNE_VERSION_LT(DUNE_LOCALFUNCTIONS,2,10)
      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);
#endif
      FieldMatrix<C, N, N> vandermonde;
      FieldVector<C, N> fEvaluations;
      int n=0;
      auto vandermondeRow = vandermonde.begin();
      for (int j=0; j<=k; j++)
        for (int i=0; i<=k-j; i++)
        {
          x[0] = static_cast<D>(i)/static_cast<D>(kdiv);
          x[1] = static_cast<D>(j)/static_cast<D>(kdiv);
          fEvaluations[n] = f(x);
          detail::Pk2DBernsteinBasis<D, R, k>
              ::fillVectorOfEvaluations(x, vandermondeRow->begin());
          ++n;
          ++vandermondeRow;
        }
      assert(vandermondeRow == vandermonde.end());
      // TODO: Solving the Vandermonde system could be sped up by
      //       precomputing the LU decomposition
      FieldVector<C, N> interpolation;
      vandermonde.solve(interpolation, fEvaluations);
      out.resize(N);
      using std::cbegin;
      using std::cend;
      std::copy(cbegin(interpolation), cend(interpolation), begin(out));
    }

  };
}

#endif
