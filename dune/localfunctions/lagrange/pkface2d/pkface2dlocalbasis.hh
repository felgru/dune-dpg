// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PKFACE2DLOCALBASIS_HH
#define DUNE_PKFACE2DLOCALBASIS_HH

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Lagrange shape functions of arbitrary order on the reference triangle.

         Lagrange shape functions of arbitrary order have the property that
         \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.
         \tparam k Polynomial order.

         \nosubgrouping
   */
  template<class D, class R, unsigned int k>
  class PkFace2DLocalBasis
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

    /** \brief Export the number of degrees of freedom */
    enum {N = 3*(k+1)};
    /** \brief Export the element order
       OS: Surprising that we need to export this both statically and dynamically!
     */
    enum {O = k};

    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,2> > Traits;

    //! \brief Standard constructor
    PkFace2DLocalBasis ()
    {}

    //! \brief number of shape functions
    unsigned int size () const
    {
      return N;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(N);
      // specialization for k==0, not clear whether that is needed
      if (k==0) {
        out[0] = 1;
        return;
      }
      out.resize(size());
      for (size_t l =0; l<size(); l++)
      {
        out[l]=0;     //TODO inefficient
      }
      if (in[1]<10e-10 and in[1]>-10e-10)
      {
        for (size_t l=0; l<=k; l++)
        {
          out[l]= p(l, in[0]);
        }
      }
      if (in[0]<10e-10 and in[0]>-10e-10)
      {
        for (size_t l=0; l<=k; l++)
        {
          out[l+(k+1)]= p(l, in[1]);
        }
      }
      if ((in[0]+in[1])<1+10e-10 and (in[0]+in[1])>1-10e-10)
      {
        for (size_t l=0; l<=k; l++)
        {
          out[l+2*(k+1)]= p(l, (std::sqrt((((1-in[0])*(1-in[0]))+(in[1]*in[1]))/2)));
        }
      }
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& x,       // position
                      std::vector<typename Traits::JacobianType>& out) const                        // return value
    {
      out.resize(N);
      DUNE_THROW(Dune::NotImplemented, "Jacobian does not make sense for face-functions");
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return k;
    }
  };

}
#endif
