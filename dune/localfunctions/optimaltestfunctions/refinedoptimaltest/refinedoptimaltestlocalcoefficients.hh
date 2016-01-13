// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_REFINED_OPTIMALTEST_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_REFINED_OPTIMALTEST_LOCALCOEFFICIENTS_HH

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
   * \tparam d Dimension of the reference cube
   */
  template<int d, int level>
  class RefinedOptimalTestLocalCoefficients {

  public:
    //! \brief Default constructor
    RefinedOptimalTestLocalCoefficients ()
    {
      li=nullptr;
    }

    RefinedOptimalTestLocalCoefficients (std::vector<LocalKey>* localKeyList)
    {
      li=localKeyList;
    }

    //! number of coefficients
    std::size_t size () const
    {
      //TODO: unsauber wenn li=nullptr
      DUNE_THROW(Dune::NotImplemented, "Be Careful - local keys not properly implemented");
      return li->size();
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      DUNE_THROW(Dune::NotImplemented, "Be Careful - local keys not properly implemented");
      return (*li)[i];
    }

  private:
    std::vector<LocalKey>* li;
  };

}

#endif
