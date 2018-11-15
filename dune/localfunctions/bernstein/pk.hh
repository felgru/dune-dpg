// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* vim: set ai expandtab sw=4 ts=4: */
#ifndef DUNE_BERNSTEIN_PK_LOCALFINITEELEMENT_HH
#define DUNE_BERNSTEIN_PK_LOCALFINITEELEMENT_HH

#include "pk2d.hh"

namespace Dune
{

  /** \brief General Bernstein finite element with arbitrary dimension and polynomial order
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam d dimension of the reference element
   * \tparam k polynomial order
   */
  template<class D, class R, int d, int k>
  class BernsteinPkLocalFiniteElement
  {
  public:
    BernsteinPkLocalFiniteElement() = default;

    /** Constructor for variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...k+1
     */
    BernsteinPkLocalFiniteElement(const unsigned int vertexmap[k+1])
    {}
  };

  /** \brief General Bernstein finite element -- specialization for a 2d reference element
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam k polynomial order
   */
  template<class D, class R, int k>
  class BernsteinPkLocalFiniteElement<D, R, 2, k>
    : public BernsteinPk2DLocalFiniteElement<D, R, k>
  {
  public:
    BernsteinPkLocalFiniteElement() = default;

    BernsteinPkLocalFiniteElement(const unsigned int vertexmap[3]) :
      BernsteinPk2DLocalFiniteElement<D, R, k>(vertexmap)
    {}
  };

}

#endif
