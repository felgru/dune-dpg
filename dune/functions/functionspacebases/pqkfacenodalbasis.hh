// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKFACENODALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKFACENODALBASIS_HH

#include <dune/common/deprecated.hh>
#include <dune/functions/functionspacebases/lagrangefacebasis.hh>


namespace Dune {
namespace Functions {



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   PQkFacePreBasis
//   PQkFaceNodeIndexSet
//   PQkFaceNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename R=double>
using PQkFaceNode = LagrangeFaceNode<GV, k, R>;

template<typename GV, int k, class MI, typename R=double>
using PQkFaceNodeIndexSet = LagrangeFaceNodeIndexSet<GV, k, MI, R>;

template<typename GV, int k, class MI, typename R=double>
using PQkFacePreBasis = LagrangeFacePreBasis<GV, k, MI, R>;



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar k-th-order Lagrangean finite element space
 *         defined on the faces of a cell
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k, typename R=double>
using PQkFaceNodalBasis
DUNE_DEPRECATED_MSG("PQkFaceNodalBasis is deprecated. Use LagrangeFaceBasis from <dune/functions/functionspacebases/lagrangefacebasis.hh> instead.")
    = LagrangeFaceBasis<GV, k, R>;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKFACENODALBASIS_HH
