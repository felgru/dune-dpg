// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_REFINED_OPTIMALTESTLOCALINTERPOLATION_HH
#define DUNE_LOCALFUNCTIONS_REFINED_OPTIMALTESTLOCALINTERPOLATION_HH

#include <dune/common/fvector.hh>
#include <dune/common/power.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>


namespace Dune
{
  /** \todo Please doc me! */
  template<int d, int level, class LB>
  class RefinedOptimalTestLocalInterpolation
  {

  public:

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      DUNE_THROW(Dune::NotImplemented, "No interpolation method for RefinedOptimalTestLocalFiniteElement available yet!");
    }
  };

}


#endif
