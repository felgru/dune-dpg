// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKDGREFINEDDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKDGREFINEDDGBASIS_HH

#include <dune/common/deprecated.hh>
#include <dune/functions/functionspacebases/lagrangedgrefineddgbasis.hh>


namespace Dune {
namespace Functions {


// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   PQkDGRefinedDGPreBasis
//   PQkDGRefinedDGNodeIndexSet
//   PQkDGRefinedDGNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
template<typename GV, int level, int k>
using PQkDGRefinedDGNode = LagrangeDGRefinedDGNode<GV, level, k>;

template<typename GV, int level, int k, class MI>
using PQkDGRefinedDGNodeIndexSet
    = LagrangeDGRefinedDGNodeIndexSet<GV, level, k, MI>;
#else
template<typename GV, int level, int k, typename TP>
using PQkDGRefinedDGNode = LagrangeDGRefinedDGNode<GV, level, k, TP>;

template<typename GV, int level, int k, class MI, class TP>
using PQkDGRefinedDGNodeIndexSet
    = LagrangeDGRefinedDGNodeIndexSet<GV, level, k, MI, TP>;
#endif

// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis of a scalar k-th-order Lagrangean-DG finite element space
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int level, int k>
using PQkDGRefinedDGBasis
DUNE_DEPRECATED_MSG("PQkDGRefinedDGBasis is deprecated. Use LagrangeDGRefinedDGBasis from <dune/functions/functionspacebases/lagrangedgrefineddgbasis.hh> instead.")
    = LagrangeDGRefinedDGBasis<GV, level, k>;


} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKDGREFINEDDGBASIS_HH
