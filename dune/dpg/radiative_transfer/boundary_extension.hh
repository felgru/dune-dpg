// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
//
// The code in Dune::detail::harmonic_extension has been adapted from
// examples/poisson-pq2.cc from dune-functions. Its original authors
// are Carsten Gräser, Christian Engwer, Oliver Sander and Steffen Müthing.
//
// The dune-functions library, headers and test programs are copyrighted
// free software. You can use, modify and/or redistribute it under the terms
// of either one of the two following licenses:
//
// * The GNU Lesser General Public License as published by the Free Software
//   Foundation, either Version 3 of the license or (at your option) any later
//   version. You can find a copy of the GNU Lesser General Public License,
//   Version 3, in the files GPL-3 and LGPL-3 or at
//   <http://www.gnu.org/licenses/lgpl-3.0>.
//
// * Version 2 of the GNU General Public License as published by the Free
//   Software Foundation, with the following special exception for linking
//   and compiling against the dune-functions library, the so-called
//   "runtime exception":
//
//     As a special exception, you may use the dune-functions source files as
//     part of a software library or application without restriction.
//     Specifically, if other files instantiate templates or use macros or
//     inline functions from one or more of the dune-functions source files,
//     or you compile one or more of the dune-functions source files and link
//     them with other files to produce an executable, this does not by
//     itself cause the resulting executable to be covered by the GNU General
//     Public License.  This exception does not however invalidate any
//     other reasons why the executable file might be covered by the GNU
//     General Public License.
//
//   This license is intended to be similar to the GNU Lesser General Public
//   License, Version 2, which by itself isn't suitable for a template library.
//   You can find a copy of the GNU General Public License, Version 2, in the
//   file GPL-2 or at <http://www.gnu.org/licenses/gpl-2.0>.

#ifndef DUNE_DPG_RADIATIVE_TRANSFER_BOUNDARY_EXTENSION_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_BOUNDARY_EXTENSION_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/functions/localindexsetiteration.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>

namespace Dune {
namespace detail {
namespace harmonic_extension {

// Compute the stiffness matrix for a single element
template <class LocalView, class Matrix>
void getLocalMatrix(const LocalView& localView, Matrix& elementMatrix)
{
  using Element = typename LocalView::Element;
  const Element& element = localView.element();

  constexpr int dim = Element::dimension;
  const auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& localFiniteElement = localView.tree().finiteElement();

  elementMatrix.setSize(localFiniteElement.size(),
                        localFiniteElement.size());
  elementMatrix = 0;      // fills the entire matrix with zeroes

  // TODO: This quadrature order seems way too large.
  const int order = 2*(dim*localFiniteElement.localBasis().order()-1);
  const QuadratureRule<double, dim>& quad
      = QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (size_t pt=0; pt < quad.size(); pt++) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The transposed inverse Jacobian of the map from the reference
    // element to the element
    const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPos);

    // The gradients of the shape functions on the reference element
    std::vector<FieldMatrix<double,1,dim> > referenceGradients;
    localFiniteElement.localBasis().evaluateJacobian(quadPos,
                                                     referenceGradients);

    // Compute the shape function gradients on the real element
    std::vector<FieldVector<double,dim> > gradients(referenceGradients.size());
    for (size_t i=0; i<gradients.size(); i++)
      jacobian.mv(referenceGradients[i][0], gradients[i]);

    // Compute the actual matrix entries
    for (size_t i=0; i<elementMatrix.N(); i++)
      for (size_t j=0; j<elementMatrix.M(); j++)
        elementMatrix[i][j] += ( gradients[i] * gradients[j] )
                              * quad[pt].weight() * integrationElement;

  }
}

// Get the occupation pattern of the stiffness matrix
template <class FEBasis>
void getOccupationPattern(const FEBasis& feBasis, MatrixIndexSet& nb)
{
  // Total number of feBasis' DoFs
  const auto n = feBasis.size();

  nb.resize(n, n);

  // A view on the FE basis on a single element
  auto localView = feBasis.localView();
  auto localIndexSet = feBasis.localIndexSet();

  // Loop over all leaf elements
  const auto gridView = feBasis.gridView();
  for(const auto& e : elements(gridView))
  {
    localView.bind(e);
    localIndexSet.bind(localView);

    using MultiIndex = typename decltype(localIndexSet)::MultiIndex;

    auto fillOccupationPattern = [&](size_t /* i */, MultiIndex gi)
        {
          auto fillOccupationPatternInner
            = [&](size_t /* j */, MultiIndex gj)
              {
                // Add a nonzero entry to the matrix
                nb.add(gi[0],
                       gj[0]);
              };
          iterateOverLocalIndexSet(
              localIndexSet,
              fillOccupationPatternInner,
              [](size_t j){},
              [&](size_t j, MultiIndex gj, double /* wj */)
              {
                fillOccupationPatternInner(j, gj);
              }
          );
        };
    iterateOverLocalIndexSet(
        localIndexSet,
        fillOccupationPattern,
        [](size_t i){},
        [&](size_t i, MultiIndex gi, double /* wi */)
        {
          fillOccupationPattern(i, gi);
        }
    );
  }
}


/** \brief Assemble the Laplace stiffness matrix on the given grid view */
template <class FEBasis, class Function>
void assembleLaplaceSystem(const FEBasis& feBasis,
                           BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                           BlockVector<FieldVector<double,1> >& rhs,
                           Function&& dirichletValueFunction)
{
  const auto gridView = feBasis.gridView();

  {
    // MatrixIndexSets store the occupation pattern of a sparse matrix.
    // They are not particularly efficient, but simple to use.
    MatrixIndexSet occupationPattern;
    getOccupationPattern(feBasis, occupationPattern);

    occupationPattern.exportIdx(matrix);
  }

  // Set all entries to zero
  matrix = 0;
  rhs = 0;

  // A view on the FE basis on a single element
  auto localView = feBasis.localView();
  auto localIndexSet = feBasis.localIndexSet();

  // A loop over all elements of the grid
  for(const auto& e : elements(gridView))
  {
    // Bind the local FE basis view to the current element
    localView.bind(e);
    localIndexSet.bind(localView);

    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    getLocalMatrix(localView, elementMatrix);

    // Add element stiffness matrix onto the global stiffness matrix
    addToGlobalMatrix(
        localIndexSet,
        localIndexSet,
        [&elementMatrix](size_t i, size_t j) -> auto {
          return elementMatrix[i][j];
        },
        [&](auto gi, auto gj) -> auto& {
          return matrix[gi[0]][gj[0]];
        }
    );
  }

  ////////////////////////////////////////////
  //   Modify Dirichlet rows
  ////////////////////////////////////////////

  // Determine Dirichlet dofs
  std::vector<bool> dirichletNodes;
  // Contribution of the boundary to the rhs
  std::vector<double> rhsContrib;

  BoundaryTools::getBoundaryMask(feBasis, dirichletNodes);

  BoundaryTools::getBoundaryValue(feBasis,
                                  rhsContrib,
                                  dirichletValueFunction);

  // loop over the matrix rows
  for (size_t i=0; i<matrix.N(); i++) {
    if (dirichletNodes[i]) {
      auto cIt          = matrix[i].begin();
      const auto cEndIt = matrix[i].end();
      // loop over nonzero matrix entries in current row
      for (; cIt!=cEndIt; ++cIt)
        *cIt = (i==cIt.index()) ? 1. : 0.;

      rhs[i] = rhsContrib[i];
    }
  }
}
}} // end namespace detail::harmonic_extension

/**
 *
 * \param dirichletValueFunction  a boundary value function in x
 */
template<class Function, class GlobalBasis>
BlockVector<FieldVector<double,1>> harmonic_extension_of_boundary_values(
    const Function& dirichletValueFunction,
    const GlobalBasis& feBasis)
{
  using Vector = BlockVector<FieldVector<double,1>>;
  using Matrix = BCRSMatrix<FieldMatrix<double,1,1>>;

  Matrix stiffnessMatrix;
  Vector rhs(feBasis.size());

  Dune::detail::harmonic_extension::
      assembleLaplaceSystem(feBasis, stiffnessMatrix, rhs,
                            dirichletValueFunction);

  ////////////////////////////
  //   Compute solution
  ////////////////////////////

  // Technicality:  turn the matrix into a linear operator
  MatrixAdapter<Matrix,Vector,Vector> op(stiffnessMatrix);

  // Sequential incomplete LU decomposition as the preconditioner
#if DUNE_VERSION_NEWER(DUNE_ISTL,2,6)
  SeqILU<Matrix,Vector,Vector> ilu0(stiffnessMatrix,1.0);
#else
  SeqILU0<Matrix,Vector,Vector> ilu0(stiffnessMatrix,1.0);
#endif

  Vector extension(feBasis.size());
  extension = 0.;

  // Preconditioned conjugate-gradient solver
  CGSolver<Vector> cg(op,
                      ilu0, // preconditioner
                      1e-4, // desired residual reduction factor
                      50,   // maximum number of iterations
                      0);   // verbosity of the solver

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  cg.apply(extension, rhs, statistics);

  return extension;
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_BOUNDARY_EXTENSION_HH
