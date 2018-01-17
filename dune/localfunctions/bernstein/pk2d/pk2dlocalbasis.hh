// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BERNSTEIN_PK2DLOCALBASIS_HH
#define DUNE_BERNSTEIN_PK2DLOCALBASIS_HH

#include <functional>
#include <numeric>

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

#include "bernsteinbasisevaluation.hh"

namespace Dune
{
  /**@ingroup LocalBasisImplementation
         \brief Bernstein shape functions of arbitrary order on the reference triangle.

         Bernstein shape functions of arbitrary order have the property that
         \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.
         \tparam k Polynomial order.

         \nosubgrouping
   */
  template<class D, class R, unsigned int k>
  class BernsteinPk2DLocalBasis
  {
    typedef Dune::FieldVector<D,3> BarycentricCoordinate;
  public:

    /** \brief Export the number of degrees of freedom */
    enum {N = (k+1)*(k+2)/2};

    /** \brief Export the element order
       OS: Surprising that we need to export this both statically and dynamically!
     */
    enum {O = k};

    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,2> > Traits;

    //! \brief Standard constructor
    BernsteinPk2DLocalBasis () = default;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return N;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& x,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(N);
      detail::Pk2DBernsteinBasis<D, R, k>
          ::fillVectorOfEvaluations(x, out.begin());
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& x,
                      std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(N);
      for(auto& o : out) {
        o[0][0] = 0; o[0][1] = 0;
      }

      // specialization for k==0, not clear whether that is needed
      if (k==0) {
        return;
      }

      std::vector<typename Traits::RangeType> bernsteinPolynomials(k*(k+1)/2);
      detail::Pk2DBernsteinBasis<D, R, k-1>
          ::fillVectorOfEvaluations(x, bernsteinPolynomials.begin());
      std::for_each(bernsteinPolynomials.begin(), bernsteinPolynomials.end(),
          [](auto& b) { b *= k; });

      const FieldVector<R, 3> v0{1, 0, -1};
      const FieldVector<R, 3> v1{0, 1, -1};

      // iterating over Bernstein polynomials of degree k-1 here
      int n=0;
      for (unsigned int i0=0; i0<k; i0++)
      {
        for (unsigned int i1=0; i1<k-i0; i1++)
        {
          const unsigned int i2 = k-1-i0-i1;

          for (unsigned char coordinate=0; coordinate<3; coordinate++) {
            std::array<unsigned int, 3> outIndex{i0, i1, i2};
            outIndex[coordinate] += 1;
            const unsigned int nOut = index(outIndex);

            // x_0 derivative
            out[nOut][0][0] += v0[coordinate] * bernsteinPolynomials[n];

            // x_1 derivative
            out[nOut][0][1] += v1[coordinate] * bernsteinPolynomials[n];
          }

          n++;
        }
      }
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

      switch (totalOrder)
      {
        case 0:
          evaluateFunction(in,out);
          break;
        default:
          DUNE_THROW(NotImplemented,
                     "Desired derivative order is not implemented.");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return k;
    }

  private:
    static inline unsigned int index(std::array<unsigned int, 3> iBarycentric)
    {
      const unsigned int i = iBarycentric[0];
      const unsigned int j = iBarycentric[1];
      return (k+1)*(k+2)/2 - (k-i+1)*(k-i+2)/2 + j;
    }
  };

}
#endif
