// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_PKDGSUBSAMPLEDLOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_PKDGSUBSAMPLEDLOCALBASIS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/power.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/lagrange/pk2d/pk2dlocalbasis.hh>


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

    using SubBasis = Pk2DLocalBasis<D, R, k>;

  public:
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,
        Dune::FieldVector<R,1>, Dune::FieldMatrix<R,1,2>, 2> Traits;

    //! \brief number of shape functions
    constexpr unsigned int size () const
    {
      return n;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      for(size_t i=0, i_max=size(); i<i_max; i++)
        out[i] = 0;
      // In which section of the subgrid are we?
      /* TODO: This does not work if we end up exactly on a vertex or an edge. */
      Dune::FieldVector<int,d> section;
      size_t sectionOffset = 0;
      typename Traits::DomainType inInSection;
      size_t line = (size_t)(in[1]*s);
      sectionOffset = (size_t)(in[0]*s);
      inInSection[0] = in[0]*s - sectionOffset;
      inInSection[1] = in[1]*s - line;
      sectionOffset *= 2;
      sectionOffset += (2*s-line)*line; // = s^2 - (s-line)^2
      bool mirrored = inInSection[0] + inInSection[1] > 1;
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
     * \param in position where to evaluate
     * \param out The return value
     */
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,
                      std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());
      for(size_t i=0, i_max=size(); i<i_max; i++)
        out[i] = 0;
      // In which section of the subgrid are we?
      /* TODO: This does not work if we end up exactly on a vertex or an edge. */
      Dune::FieldVector<int,d> section;
      size_t sectionOffset = 0;
      typename Traits::DomainType inInSection;
      size_t line = (size_t)(in[1]*s);
      sectionOffset = (size_t)(in[0]*s);
      inInSection[0] = in[0]*s - sectionOffset;
      inInSection[1] = in[1]*s - line;
      sectionOffset *= 2;
      sectionOffset += (2*s-line)*line; // = s^2 - (s-line)^2
      bool mirrored = inInSection[0] + inInSection[1] > 1;
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
            out[sectionOffset+i][0][d] = -s*outInSection[i][0][1-d];
      } else {
        for(size_t i=0, i_max=outInSection.size(); i<i_max; ++i)
          for(unsigned int d=0; d<2; ++d)
            out[sectionOffset+i][0][d] = s*outInSection[i][0][d];
      }
    }

    /** \brief Evaluate derivative in a given direction
     * \param [in]  direction The direction to derive in
     * \param [in]  in        Position where to evaluate
     * \param [out] out       The return value
     */
    template<int diffOrder>
    inline void evaluate(
      const std::array<int,diffOrder>& direction,
      const typename Traits::DomainType& in,
      std::vector<typename Traits::RangeType>& out) const
    {
      DUNE_THROW(Dune::NotImplemented, "Evaluation of arbitrary "
              "derivatives of PkDGSubsampled2DLocalBasis "
              "not implemented, yet.");
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return k;
    }
  };
}

#endif
