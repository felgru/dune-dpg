// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QKDISCONTUNUOUS_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_QKISCONTUNUOUS_LOCALCOEFFICIENTS_HH

#include <array>
#include <cassert>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/power.hh>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{
  /** \brief Attaches a shape function to an entity
   *
   * \tparam k Polynomial order
   * \tparam d Dimension of the reference cube
   */
  template<int k, int d>
  class QkDiscontinuousLocalCoefficients {

  public:
    //! \brief Default constructor
    QkDiscontinuousLocalCoefficients () : li(StaticPower<k+1,d>::power)
    {
      // Set up array of codimension-per-dof-number,
      // index vector (the index of the dof in the set of dofs of a given subentity) and
      // sub entity number
      std::vector<unsigned int> codim(li.size());
      std::vector<unsigned int> index(size());
      std::vector<unsigned int> subEntity(li.size());

      for (std::size_t i=0; i<codim.size(); i++)
      {
        codim[i] = 0;
        index[i] = i;
        subEntity[i] = 0;
      }

      for (size_t i=0; i<li.size(); i++)
        li[i] = LocalKey(subEntity[i], codim[i], index[i]);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return StaticPower<k+1,d>::power;
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
