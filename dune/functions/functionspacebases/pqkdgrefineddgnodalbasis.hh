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

template<typename GV, int level, int k, typename R=double>
using PQkDGRefinedDGNode = LagrangeDGRefinedDGNode<GV, level, k, R>;

template<typename GV, int level, int k, class MI, typename R=double>
using PQkDGRefinedDGNodeIndexSet
    = LagrangeDGRefinedDGNodeIndexSet<GV, level, k, MI, R>;

// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis of a scalar k-th-order Lagrangean-DG finite element space
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int level, int k, typename R=double>
using PQkDGRefinedDGBasis
DUNE_DEPRECATED_MSG("PQkDGRefinedDGBasis is deprecated. Use LagrangeDGRefinedDGBasis from <dune/functions/functionspacebases/lagrangedgrefineddgbasis.hh> instead.")
    = LagrangeDGRefinedDGBasis<GV, level, k, R>;


} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKDGREFINEDDGBASIS_HH
