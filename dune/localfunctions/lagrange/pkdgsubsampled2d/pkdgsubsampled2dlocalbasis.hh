// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_PKDGSUBSAMPLEDLOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_PKDGSUBSAMPLEDLOCALBASIS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/power.hh>
#include <dune/common/version.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#if DUNE_VERSION_GTE(DUNE_LOCALFUNCTIONS,2,7)
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#else
#include <dune/localfunctions/lagrange/pk2d/pk2dlocalbasis.hh>
#endif

#include <numeric>


namespace Dune
{
  /**@ingroup LocalBasisImplementation
   * \brief Lagrange shape functions of order k on the reference triangle.
   *
   * Also known as \f$P^k\f$.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   * \tparam s number of samples
   * \tparam k Polynomial degree
   *
   * \nosubgrouping
   */
  template<class D, class R, unsigned int s, unsigned int k>
  class PkDGSubsampled2DLocalBasis
  {
    enum { n = (k+1)*(k+2)/2*StaticPower<s,2>::power };
    enum { d = 2 };

#if DUNE_VERSION_GTE(DUNE_LOCALFUNCTIONS,2,7)
    using SubBasis = typename LagrangeSimplexLocalFiniteElement<D, R, 2, k>
                                                    ::Traits::LocalBasisType;
#else
    using SubBasis = Pk2DLocalBasis<D, R, k>;
#endif

  public:
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,
        Dune::FieldVector<R,1>, Dune::FieldMatrix<R,1,2> > Traits;

    //! \brief number of shape functions
    constexpr unsigned int size () const
    {
      return n;
    }

    /** \brief Evaluate all shape functions
     *
     * \note This might not work as expected if you evaluate exactly on
     *       the edge between two subcells.
     */
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      for(size_t i=0, i_max=size(); i<i_max; i++)
        out[i] = 0;
      // In which section of the subgrid are we?
      typename Traits::DomainType inInSection;
      size_t line          = static_cast<size_t>(in[1]*s);
      if(line>=s) line = s-1;
      size_t sectionOffset = static_cast<size_t>(in[0]*s);
      if(sectionOffset>=s-line) sectionOffset = s-line-1;
      inInSection[0] = in[0]*s - sectionOffset;
      inInSection[1] = in[1]*s - line;
      const bool mirrored = (inInSection[0] + inInSection[1] > 1)
                          && (sectionOffset < s-line-1);
      sectionOffset *= 2;
      sectionOffset += (2*s-line)*line; // = s^2 - (s-line)^2
      if(mirrored) {
        sectionOffset += 1;
        auto tmp = inInSection[0];
        inInSection[0] = 1 - inInSection[1];
        inInSection[1] = 1 - tmp;
      }
      sectionOffset *= (k+1)*(k+2)/2;

      std::vector<typename Traits::RangeType> outInSection;
      SubBasis().evaluateFunction(inInSection, outInSection);
      for(size_t i=0, i_max=outInSection.size(); i<i_max; ++i)
        out[sectionOffset+i] = outInSection[i];
    }

    /** \brief Evaluate Jacobian of all shape functions
     *
     * \param in position where to evaluate
     * \param out The return value
     *
     * \note This might not work as expected if you evaluate exactly on
     *       the edge between two subcells.
     */
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,
                      std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());
      for(size_t i=0, i_max=size(); i<i_max; i++)
        out[i] = 0;
      // In which section of the subgrid are we?
      typename Traits::DomainType inInSection;
      size_t line          = static_cast<size_t>(in[1]*s);
      if(line>=s) line = s-1;
      size_t sectionOffset = static_cast<size_t>(in[0]*s);
      if(sectionOffset>=s-line) sectionOffset = s-line-1;
      inInSection[0] = in[0]*s - sectionOffset;
      inInSection[1] = in[1]*s - line;
      const bool mirrored = (inInSection[0] + inInSection[1] > 1)
                          && (sectionOffset < s-line-1);
      sectionOffset *= 2;
      sectionOffset += (2*s-line)*line; // = s^2 - (s-line)^2
      if(mirrored) {
        sectionOffset += 1;
        size_t tmp = inInSection[0];
        inInSection[0] = 1 - inInSection[1];
        inInSection[1] = 1 - tmp;
      }
      sectionOffset *= (k+1)*(k+2)/2;

      std::vector<typename Traits::JacobianType> outInSection;
      SubBasis().evaluateJacobian(inInSection, outInSection);
      // Transform Jacobian if mirrored
      if(mirrored) {
        for(size_t i=0, i_max=outInSection.size(); i<i_max; ++i)
          for(unsigned int d=0; d<2; ++d)
            out[sectionOffset+i][0][d]
              = -static_cast<double>(s)*outInSection[i][0][1-d];
      } else {
        for(size_t i=0, i_max=outInSection.size(); i<i_max; ++i)
          for(unsigned int d=0; d<2; ++d)
            out[sectionOffset+i][0][d] = s*outInSection[i][0][d];
      }
    }

    inline void partial(const std::array<unsigned int,d>& order,
                        const typename Traits::DomainType& in,
                        std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else {
        DUNE_THROW(Dune::NotImplemented, "Evaluation of arbitrary "
              "derivatives of PkDGSubsampled2DLocalBasis "
              "not implemented, yet.");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return k;
    }
  };
}

#endif
