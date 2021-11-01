// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QKFACE2D_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_QKFACE2D_LOCALCOEFFICIENTS_HH

#include <array>
#include <cassert>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{
  /** \brief Attaches a shape function to an entity
   *
   * \tparam k Polynomial order
   * \tparam d Dimension of the reference cube
   */
  template<int k>
  class QkFace2DLocalCoefficients {

      // LocalKey: entity number , entity codim, dof indices within each entity
      /* edge and vertex numbering
       2----3----3
       |         |
       |         |
       0         1
       |         |
       |         |
       0----2----1
       */


  public:
    //! \brief Default constructor
    QkFace2DLocalCoefficients () : li(4*(k+1))
    {
      // codimension-per-dof-number = 1 for all indices
      for (size_t l = 0; l<=k; l++)
      {
        li[l] = LocalKey(0, 1, l);
        li[l+k+1] = LocalKey(1, 1, l);
        li[l+2*(k+1)] = LocalKey(2, 1, l);
        li[l+3*(k+1)] = LocalKey(3, 1, l);
      }

    }

    //! number of coefficients
    std::size_t size () const
    {
      return 4*(k+1);
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;
  };

}

#endif
