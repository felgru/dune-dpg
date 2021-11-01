// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QKFACE2DLOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_QKFACE2DLOCALBASIS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>

#include <numeric>


namespace Dune
{
  /**@ingroup LocalBasisImplementation
     \brief Lagrange shape functions of order k on the reference cube.

     Also known as \f$Q^k\f$.

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.
     \tparam k Polynomial degree

     \nosubgrouping
   */
  template<class D, class R, int k>
  class QkFace2DLocalBasis
  {
    // ith Lagrange polynomial of degree k in one dimension
    static R p (int i, D x)
    {
      R result(1.0);
      for (int j=0; j<=k; j++)
        if (j!=i) result *= (k*x-j)/(i-j);
      return result;
    }

  public:
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,Dune::FieldMatrix<R,1,2> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return (4*(k+1));
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      for (size_t l =0; l<size(); l++)
      {
        out[l]=0;     //TODO inefficient
      }
      if (in[0]<10e-10 and in[0]>-10e-10)
      {
        for (size_t l=0; l<=k; l++)
        {
          out[l]= p(l, in[1]);
        }
      }
      else if (in[0]<1+10e-10 and in[0]>1-10e-10)
      {
        for (size_t l=0; l<=k; l++)
        {
          out[l+k+1]= p(l, in[1]);
        }
      }
      if (in[1]<10e-10 and in[1]>-10e-10)
      {
        for (size_t l=0; l<=k; l++)
        {
          out[l+2*(k+1)]= p(l, in[0]);
        }
      }
      else if (in[1]<1+10e-10 and in[1]>1-10e-10)
      {
        for (size_t l=0; l<=k; l++)
        {
          out[l+3*(k+1)]= p(l, in[0]);
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
      out.resize(size());
      DUNE_THROW(Dune::NotImplemented, "Jacobian does not make sense for face-functions");
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    void partial(const std::array<unsigned int,2>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else {
        DUNE_THROW(Dune::NotImplemented,
            "partial only implemented for derivatives of order 0!");
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
