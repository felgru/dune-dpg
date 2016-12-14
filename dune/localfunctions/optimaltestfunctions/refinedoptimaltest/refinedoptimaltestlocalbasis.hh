// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_REFINED_OPTIMALTESTLOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_REFINED_OPTIMALTESTLOCALBASIS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/power.hh>

#include <dune/geometry/type.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>

#include <numeric>


namespace Dune
{
  /**@ingroup LocalBasisImplementation
     \brief Optimal test space shape functions on the reference element.

     The shape functions are defined by the shape functions of the given
     test search space and a coefficient matrix.

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.
     \tparam d Dimension of the cube
     \tparam level Level of local refinement
     \tparam TestSearchLocalBasis

     \nosubgrouping
   */
  template<class D, class R, int d,
           int level, class RefinementConstants, class TestSearchLocalBasis>
  class RefinedOptimalTestLocalBasis
  {
    typedef Matrix<FieldMatrix<double,1,1> > MatrixType;

  public:
    typedef LocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,Dune::FieldMatrix<R,1,d> > Traits;

    RefinedOptimalTestLocalBasis
        (const TestSearchLocalBasis& testSearchLocalBasis,
         const GeometryType& gt, MatrixType& coeffMat, size_t offset)
        : testSearchLocalBasis(testSearchLocalBasis),
          coefficientMatrix(coeffMat),
          offset(offset),
          numTestSearchDOFs(testSearchLocalBasis.size())
      {
        assert(coeffMat.N() >= offset + numTestSearchDOFs);

        unsigned int order = testSearchLocalBasis.order();
        if(gt.isTriangle())
          subElementSize = (order+1)*(order+2)/2;
        else if(gt.isQuadrilateral())
          subElementSize = (order+1)*(order+1);
        else
          DUNE_THROW(Dune::NotImplemented,
            "RefinedOptimalTestLocalBasis not implemented for GeometryType "
            << gt.id() << "!");
      }

    //! \brief number of shape functions
    unsigned int size () const
    {
      return coefficientMatrix.M();
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (unsigned int subElement,
                                  const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      std::vector<FieldVector<double,1> > valuesTestSearchLocalBasis(numTestSearchDOFs);
      testSearchLocalBasis.evaluateFunction(in, valuesTestSearchLocalBasis);
      size_t subElementOffset = subElement * subElementSize;
      for (size_t i=0, iMax=size(); i<iMax; i++)
      {
        out[i] = 0;
        for (size_t j=0; j<numTestSearchDOFs; ++j)
        {
          out[i] += coefficientMatrix[offset+subElementOffset+j][i][0]
                    * valuesTestSearchLocalBasis[j][0];
        }
      }
    }

    /** Evaluate Jacobian of all shape functions
     *
     * \param subElement index of the subelement on which we take the
     *                   shape functions
     * \param in position where to evaluate
     * \param out The return value
     */
    inline void
    evaluateJacobian (unsigned int subElement,
                      const typename Traits::DomainType& in,
                      std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());
      std::vector<typename Traits::JacobianType> JacobianTestSearchLocalBasis(numTestSearchDOFs);
      // TODO: Jacobian needs scaling due to refinement
      testSearchLocalBasis.evaluateJacobian(in, JacobianTestSearchLocalBasis);
      size_t subElementOffset = subElement * subElementSize;
      // Loop over all shape functions
      for (size_t i=0, iMax=size(); i<iMax; i++)
      {
        // Loop over all coordinate directions
        for (int b=0; b<d; b++)
        {
          // Initialize: the overall expression is a product
          // if j-th bit of i is set to -1, else 1
          out[i][0][b] = 0;
          for (size_t j=0; j<numTestSearchDOFs; ++j)
          {
            out[i][0][b] += coefficientMatrix[offset+subElementOffset+j][i][0]
                            * JacobianTestSearchLocalBasis[j][0][b];
          }
        }
      }
    }

    /** Evaluate partial derivatives of any order of all shape functions
     *
     * \param subElement index of the subelement on which we take the
     *                   shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    void partial(unsigned int subElement,
                 const std::array<unsigned int,d>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(subElement, in, out);
      } else {
        DUNE_THROW(Dune::NotImplemented,
            "partial only implemented for derivatives of order 0!");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
      return testSearchLocalBasis.order();
    }
  private:
    const TestSearchLocalBasis& testSearchLocalBasis;
    MatrixType& coefficientMatrix;
    size_t offset;
    size_t numTestSearchDOFs;
    size_t subElementSize;
  };
}

#endif
