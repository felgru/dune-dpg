// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QKSUBSAMPLEDLOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_QKSUBSAMPLEDLOCALBASIS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/power.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>


namespace Dune
{
  /**@ingroup LocalBasisImplementation
   * \brief Lagrange shape functions of order k on the reference cube.
   *
   * Also known as \f$Q^k\f$.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   * \tparam s number of samples
   * \tparam k Polynomial degree
   * \tparam d Dimension of the cube
   *
   * \nosubgrouping
   */
  template<class D, class R, int s, int k, int d>
  class QkSubsampledLocalBasis
  {
    enum { n = StaticPower<(s*k)+1,d>::power };

    // ith Lagrange polynomial of degree k in one dimension
    static R p (int i, D x)
    {
      R result(1.0);
      for (int j=0; j<=k; j++)
        if (j!=i) result *= (k*x-j)/(i-j);
      return result;
    }

    // derivative of ith Lagrange polynomial of degree k in one dimension
    static R dp (int i, D x)
    {
      R result(0.0);

      for (int j=0; j<=k; j++)
        if (j!=i)
        {
          R prod( (k*1.0)/(i-j) );
          for (int l=0; l<=k; l++)
            if (l!=i && l!=j)
              prod *= (k*x-l)/(i-l);
          result += prod;
        }
      // We need to scale by s to counter the scaling of our lagrange
      // basis to 1/s of the unit interval.
      return s*result;
    }

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
    typedef LocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,Dune::FieldMatrix<R,1,d>, 1> Traits;

    //! \brief number of shape functions
    constexpr unsigned int size () const
    {
      return n;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      static_assert(d==2, "QKSubsampledLocalBasis only implemented in 2D!");

      out.resize(size());
      for(size_t i=0, i_max=out.size(); i<i_max; i++)
        out[i] = 0;
      // In which section of the grid are we?
      /* What if we end up exactly on a vertex or an edge? */
      Dune::FieldVector<int,d> section;
      size_t sectionIndex = 0;
      size_t sectionOffset = 0;
      typename Traits::DomainType inInSection;
      for(int i=d-1; i>=0; i--)
      {
        section[i]     = (int)(in[i]*s);
        if(section[i] >= s)
            section[i] = s-1;
        sectionIndex  *= s;
        sectionIndex  += section[i];
        sectionOffset *= (s*k+1);
        sectionOffset += section[i]*k;
        inInSection[i] = in[i]*s - section[i];
      }

      for (size_t j=0; j<k+1; j++)
      {
        size_t dofOffset = sectionOffset + j*(s*k+1);
        for (size_t i=0; i<k+1; i++)
        {
          const size_t dofIndex = dofOffset + i;
          out[dofIndex] = p(i,inInSection[0]) * p(j,inInSection[1]);
        }
      }
    }

    /** \brief Evaluate Jacobian of all shape functions
     * \param in position where to evaluate
     * \param out The return value
     */
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,
                      std::vector<typename Traits::JacobianType>& out) const
    {
      static_assert(d==2, "QKSubsampledLocalBasis only implemented in 2D!");

      out.resize(size());
      for(size_t i=0, i_max=out.size(); i<i_max; i++)
        out[i] = 0;
      // In which section of the grid are we?
      /* What if we end up exactly on a vertex or an edge? */
      Dune::FieldVector<int,d> section;
      size_t sectionIndex = 0;
      size_t sectionOffset = 0;
      typename Traits::DomainType inInSection;
      for(int i=d-1; i>=0; i--)
      {
        section[i]     = (int)(in[i]*s);
        if(section[i] >= s)
            section[i] = s-1;
        sectionIndex  *= s;
        sectionIndex  += section[i];
        sectionOffset *= (s*k+1);
        sectionOffset += section[i]*k;
        inInSection[i] = in[i]*s - section[i];
      }

      for (size_t j=0; j<k+1; j++)
      {
        size_t dofOffset = sectionOffset + j*(s*k+1);
        for (size_t i=0; i<k+1; i++)
        {
          const size_t dofIndex = dofOffset + i;

          Dune::FieldVector<int,d> multiindexInSection;
          multiindexInSection[0] = i;
          multiindexInSection[1] = j;

          // Compute Jacobian
          for (int j=0; j<d; j++)
          {
            out[dofIndex][0][j] = dp(multiindexInSection[j],
                                     inInSection[j]);

            for (int l=0; l<d; l++)
              if (l!=j)
                out[dofIndex][0][j] *= p(multiindexInSection[l],
                                         inInSection[l]);
          }
        }
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
              "derivatives of QkSubsampledLocalBasis not implemented, yet.");
      }
    }

    /** \brief Evaluate derivative in a given direction
     * \param [in]  direction The direction to derive in
     * \param [in]  in        Position where to evaluate
     * \param [out] out       The return value
     */
    template<int diffOrder>
    inline void evaluate(
      const std::array<int,1>& direction,
      const typename Traits::DomainType& in,
      std::vector<typename Traits::RangeType>& out) const
    {
      DUNE_THROW(Dune::NotImplemented, "Evaluation of arbitrary "
              "derivatives of QkSubsampledLocalBasis not implemented, yet.");
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return k;
    }
  };
}

#endif
