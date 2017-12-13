// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PKTRACE2DLOCALFINITEELEMENT_HH
#define DUNE_PKTRACE2DLOCALFINITEELEMENT_HH

#include <cstddef>

#include <dune/common/version.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include "pktrace2d/pktrace2dlocalbasis.hh"
#include "pktrace2d/pktrace2dlocalcoefficients.hh"
#include "pktrace2d/pktrace2dlocalinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, unsigned int k>
  class PkTrace2DLocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<PkTrace2DLocalBasis<D,R,k>,
        PkTrace2DLocalCoefficients<k>,
        PkTrace2DLocalInterpolation<PkTrace2DLocalBasis<D,R,k> > > Traits;

    /** \todo Please doc me !
     */
    PkTrace2DLocalFiniteElement ()
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
    { }
#else
    {
      gt.makeTriangle();
    }
#endif

    /** \todo Please doc me !
     */
    PkTrace2DLocalFiniteElement (int variant)
      : coefficients(variant)
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
    { }
#else
    {
      gt.makeTriangle();
    }
#endif

    /** Constructor for six variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...2
     */
    PkTrace2DLocalFiniteElement (const unsigned int vertexmap[3])
      : coefficients(vertexmap)
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
    { }
#else
    {
      gt.makeTriangle();
    }
#endif

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    /** \todo Please doc me !
     */
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
    static constexpr GeometryType type ()
    {
      return GeometryTypes::triangle;
    }
#else
    GeometryType type () const
    {
      return gt;
    }
#endif

  private:
    PkTrace2DLocalBasis<D,R,k> basis;
    PkTrace2DLocalCoefficients<k> coefficients;
    PkTrace2DLocalInterpolation<PkTrace2DLocalBasis<D,R,k> > interpolation;
#if not(DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6))
    GeometryType gt;
#endif
  };
#if 0
  //! Langrange finite element of arbitrary order on triangles
  /**
   * \tparam Geometry Geometry for the local to global transformation.
   * \tparam RF       Field type of the range.
   * \tparam k        Maximum polynomial order of the base functions.
   *
   * \implements FiniteElementInterface
   */
  template<class Geometry, class RF, std::size_t k>
  class PkTrace2DFiniteElement {
    typedef typename Geometry::ctype DF;
    typedef PkTrace2DLocalBasis<DF,RF,k> LocalBasis;
    typedef PkTrace2DLocalInterpolation<LocalBasis> LocalInterpolation;

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
      typedef PkTrace2DLocalCoefficients<k> Coefficients;
    };

  private:
    static const GeometryType gt;
    static const LocalBasis localBasis;
    static const LocalInterpolation localInterpolation;

    typename Traits::Basis basis_;
    typename Traits::Interpolation interpolation_;
    typename Traits::Coefficients coefficients_;

  public:
    //! construct a PkTrace2DFiniteElement
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
    PkTrace2DFiniteElement(const Geometry &geometry,
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
  PkTrace2DFiniteElement<Geometry, RF, k>::gt(GeometryType::simplex, 2);

  template<class Geometry, class RF, std::size_t k>
  const typename PkTrace2DFiniteElement<Geometry, RF, k>::LocalBasis
  PkTrace2DFiniteElement<Geometry, RF, k>::localBasis = LocalBasis();

  template<class Geometry, class RF, std::size_t k>
  const typename PkTrace2DFiniteElement<Geometry, RF, k>::LocalInterpolation
  PkTrace2DFiniteElement<Geometry, RF, k>::localInterpolation =
    LocalInterpolation();

  //! Factory for PkTrace2DFiniteElement objects
  /**
   * Constructs PkTrace2DFiniteElement objects given a geometry and a vertex
   * ordering.
   *
   * \tparam Geometry Geometry for the local to global transformation.
   * \tparam RF       Field type of the range.
   * \tparam k        Maximum polynomial order of the base functions.
   *
   * \implements FiniteElementFactoryInterface
   */
  template<class Geometry, class RF, std::size_t k>
  struct PkTrace2DFiniteElementFactory {
    typedef PkTrace2DFiniteElement<Geometry, RF, k> FiniteElement;

    //! construct PkTrace2DFiniteElementFactory
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
#endif
}

#endif
