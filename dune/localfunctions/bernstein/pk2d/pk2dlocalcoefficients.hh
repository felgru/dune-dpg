// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_BERNSTEIN_PK2DLOCALCOEFFICIENTS_HH
#define DUNE_BERNSTEIN_PK2DLOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for Bernstein Pk elements

         \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  template<unsigned int k>
  class BernsteinPk2DLocalCoefficients
  {
    enum {N = (k+1)*(k+2)/2};

  public:
    //! \brief Standard constructor
    BernsteinPk2DLocalCoefficients () : li(N)
    {
      fill_default();
    }

    //! constructor for eight variants with order on edges flipped
    BernsteinPk2DLocalCoefficients (int variant) : li(N)
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
    explicit BernsteinPk2DLocalCoefficients(const VertexMap &vertexmap) : li(N)
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
      if (k==0)
      {
        li[0] = LocalKey(0,0,0);
        return;
      }
      int n=0;
      int c=0;
      for (unsigned int i0=0; i0<=k; i0++)
      {
        for (unsigned int i1=0; i1<=k-i0; i1++)
        {
          const unsigned int i2 = k-i0-i1;
          if (i0==0 && i1==0)
          {
            li[n++] = LocalKey(0,2,0);
            continue;
          }
          if (i0==k && i1==0)
          {
            li[n++] = LocalKey(1,2,0);
            continue;
          }
          if (i0==0 && i1==k)
          {
            li[n++] = LocalKey(2,2,0);
            continue;
          }
          if (i1==0)
          {
            li[n++] = LocalKey(0,1,i0-1);
            continue;
          }
          if (i0==0)
          {
            li[n++] = LocalKey(1,1,i1-1);
            continue;
          }
          if (i2==0)
          {
            li[n++] = LocalKey(2,1,i1-1);
            continue;
          }
          li[n++] = LocalKey(0,0,c++);
        }
      }
    }
  };

}

#endif
