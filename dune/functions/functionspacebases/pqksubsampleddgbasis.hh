// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKSUBSAMPLEDDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKSUBSAMPLEDDGBASIS_HH

#include <dune/common/deprecated.hh>
#include <dune/functions/functionspacebases/lagrangesubsampleddgbasis.hh>


namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   PQkSubsampledDGPreBasis
//   PQkSubsampledDGNodeIndexSet
//   PQkSubsampledDGNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int s, int k, typename R=double>
using PQkSubsampledDGNode = LagrangeSubsampledDGNode<GV, s, k, R>;

template<typename GV, int s, int k, class MI, typename R=double>
using PQkSubsampledDGNodeIndexSet
    = LagrangeSubsampledDGNodeIndexSet<GV, s, k, MI, R>;

template<typename GV, int s, int k, class MI, typename R=double>
using PQkSubsampledDGPreBasis
= LagrangeSubsampledDGPreBasis<GV, s, k, MI, R>;


// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar k-th-order discontinuous Lagrangean finite element space
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int s, int k, typename R=double>
using PQkSubsampledDGNodalBasis
DUNE_DEPRECATED_MSG("PQkSubsampledDGNodalBasis is deprecated. Use LagrangeSubsampledDGBasis from <dune/functions/functionspacebases/lagrangesubsampleddgbasis.hh> instead.")
    = LagrangeSubsampledDGBasis<GV, s, k, R>;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKSUBSAMPLEDDGBASIS_HH
