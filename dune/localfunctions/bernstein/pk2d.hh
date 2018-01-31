// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BERNSTEIN_PK2DLOCALFINITEELEMENT_HH
#define DUNE_BERNSTEIN_PK2DLOCALFINITEELEMENT_HH

#include <cstddef>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "pk2d/pk2dlocalbasis.hh"
#include "pk2d/pk2dlocalcoefficients.hh"
#include "pk2d/pk2dlocalinterpolation.hh"

namespace Dune
{
  template<class D, class R, unsigned int k>
  class BernsteinPk2DLocalFiniteElement
  {
  public:
    typedef LocalFiniteElementTraits<BernsteinPk2DLocalBasis<D,R,k>,
        BernsteinPk2DLocalCoefficients<k>,
        BernsteinPk2DLocalInterpolation<BernsteinPk2DLocalBasis<D,R,k>>> Traits;

    BernsteinPk2DLocalFiniteElement ()
    {}

    BernsteinPk2DLocalFiniteElement (int variant) :
      coefficients(variant)
    {}

    /** Constructor for six variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...2
     */
    BernsteinPk2DLocalFiniteElement (const unsigned int vertexmap[3]) :
      coefficients(vertexmap)
    {}

    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    static constexpr GeometryType type ()
    {
      return GeometryTypes::triangle;
    }

  private:
    BernsteinPk2DLocalBasis<D,R,k> basis;
    BernsteinPk2DLocalCoefficients<k> coefficients;
    BernsteinPk2DLocalInterpolation<BernsteinPk2DLocalBasis<D,R,k> > interpolation;
  };

  //! Bernstein finite element of arbitrary order on triangles
  /**
   * \tparam Geometry Geometry for the local to global transformation.
   * \tparam RF       Field type of the range.
   * \tparam k        Maximum polynomial order of the base functions.
   *
   * \implements FiniteElementInterface
   */
  template<class Geometry, class RF, std::size_t k>
  class BernsteinPk2DFiniteElement {
    typedef typename Geometry::ctype DF;
    typedef BernsteinPk2DLocalBasis<DF,RF,k> LocalBasis;
    typedef BernsteinPk2DLocalInterpolation<LocalBasis> LocalInterpolation;

  public:
    /**
     * \implements FiniteElementInterface::Traits
     */
    struct Traits {
      typedef ScalarLocalToGlobalBasisAdaptor<LocalBasis, Geometry> Basis;
      typedef LocalToGlobalInterpolationAdaptor<
          LocalInterpolation,
          typename Basis::Traits
          > Interpolation;
      typedef BernsteinPk2DLocalCoefficients<k> Coefficients;
    };

  private:
    static const GeometryType gt;
    static const LocalBasis localBasis;
    static const LocalInterpolation localInterpolation;

    typename Traits::Basis basis_;
    typename Traits::Interpolation interpolation_;
    typename Traits::Coefficients coefficients_;

  public:
    //! construct a BernsteinPk2DFiniteElement
    /**
     * \param geometry    The geometry object to use for adaption.
     * \param vertexOrder The global ordering of the vertices within the grid,
     *                    used to determine orientation of the edges.  This
     *                    vertexOrder object must support codim=0.
     *
     * \note This class stores the reference to the geometry passed here.  Any
     *       use of this class after this references has become invalid
     *       results in undefined behaviour.  The exception is that the
     *       destructor of this class may still be called.  The information
     *       contained in the vertexOrder object is extracted and the object
     *       is no longer needed after the contructor returns.
     */
    template<class VertexOrder>
    BernsteinPk2DFiniteElement(const Geometry &geometry,
                      const VertexOrder& vertexOrder) :
      basis_(localBasis, geometry), interpolation_(localInterpolation),
      coefficients_(vertexOrder.begin(0, 0))
    { }

    const typename Traits::Basis& basis() const { return basis_; }
    const typename Traits::Interpolation& interpolation() const
    { return interpolation_; }
    const typename Traits::Coefficients& coefficients() const
    { return coefficients_; }
    const GeometryType &type() const { return gt; }
  };

  template<class Geometry, class RF, std::size_t k>
  const GeometryType
  BernsteinPk2DFiniteElement<Geometry, RF, k>::gt(GeometryTypes::simplex(2));

  template<class Geometry, class RF, std::size_t k>
  const typename BernsteinPk2DFiniteElement<Geometry, RF, k>::LocalBasis
  BernsteinPk2DFiniteElement<Geometry, RF, k>::localBasis = LocalBasis();

  template<class Geometry, class RF, std::size_t k>
  const typename BernsteinPk2DFiniteElement<Geometry, RF, k>::LocalInterpolation
  BernsteinPk2DFiniteElement<Geometry, RF, k>::localInterpolation =
    LocalInterpolation();

  //! Factory for BernsteinPk2DFiniteElement objects
  /**
   * Constructs BernsteinPk2DFiniteElement objects given a geometry and a vertex
   * ordering.
   *
   * \tparam Geometry Geometry for the local to global transformation.
   * \tparam RF       Field type of the range.
   * \tparam k        Maximum polynomial order of the base functions.
   *
   * \implements FiniteElementFactoryInterface
   */
  template<class Geometry, class RF, std::size_t k>
  struct BernsteinPk2DFiniteElementFactory {
    typedef BernsteinPk2DFiniteElement<Geometry, RF, k> FiniteElement;

    //! construct BernsteinPk2DFiniteElementFactory
    /**
     * \param geometry    The geometry object to use for adaption.
     * \param vertexOrder The global ordering of the vertices within the grid,
     *                    used to determine orientation of the edges.  This
     *                    vertexOrder object must support codim=0.
     *
     * \note The returned object stores the reference to the geometry passed
     *       here.  Any use of the returned value after this references has
     *       become invalid results in undefined behaviour.  The exception is
     *       that the destructor of this class may still be called.  The
     *       information contained in the vertexOrder object is extracted and
     *       the object is no longer needed after the contructor returns.  No
     *       reference to internal data of the factory is stored.
     */
    template<class VertexOrder>
    const FiniteElement make(const Geometry& geometry,
                             const VertexOrder& vertexOrder)
    { return FiniteElement(geometry, vertexOrder); }
  };
}

#endif
