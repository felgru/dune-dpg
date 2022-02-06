// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_PKDGSUBSAMPLED2DLOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_PKDGSUBSAMPLED2DLOCALCOEFFICIENTS_HH

#include <array>
#include <cassert>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/math.hh>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{
  /** \brief Attaches a shape function to an entity
   *
   * \tparam s number of samples
   * \tparam k Polynomial order
   */
  template<int s, int k>
  class PkDGSubsampled2DLocalCoefficients {
    enum { N = (k+1)*(k+2)/2 * power(s, 2) };

  public:
    //! \brief Standard constructor
    PkDGSubsampled2DLocalCoefficients () : li(N)
    {
      fill_default();
    }

    //! constructor for eight variants with order on edges flipped
    PkDGSubsampled2DLocalCoefficients (int variant) : li(N)
    {
      fill_default();
      bool flip[3];
      for (int i = 0; i < 3; ++i)
        flip[i] = variant & (1<<i);
      for (int i=0; i<N; i++)
        if (li[i].codim()==1 && flip[li[i].subEntity()])
          li[i].index(k-2-li[i].index());
    }

    /** Constructor for six variants with permuted vertices.

        \param vertexmap The permutation of the vertices.  This
        can for instance be generated from the global indices of
        the vertices by reducing those to the integers 0...2.  This may be any
        object which for which the expression \c vertexmap[i] is defined
        appropriately (like an array, a pointer, a std::vector, or a
        random-access iterator.
     */
    template<class VertexMap>
    explicit PkDGSubsampled2DLocalCoefficients(const VertexMap &vertexmap) : li(N)
    {
      fill_default();
      bool flip[3];
      flip[0] = vertexmap[0] > vertexmap[1];
      flip[1] = vertexmap[0] > vertexmap[2];
      flip[2] = vertexmap[1] > vertexmap[2];
      for (std::size_t i=0; i<N; i++)
        if (li[i].codim()==1 && flip[li[i].subEntity()])
          li[i].index(k-2-li[i].index());
    }

    //! number of coefficients
    std::size_t size () const
    {
      return N;
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;

    void fill_default ()
    {
      for(size_t i=0; i<N; ++i)
        li[i] = LocalKey(0,0,i);
    }
  };

}

#endif
