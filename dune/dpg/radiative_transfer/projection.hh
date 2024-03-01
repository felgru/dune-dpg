// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVATE_TRANSFER_PROJECTION_HH
#define DUNE_DPG_RADIATIVATE_TRANSFER_PROJECTION_HH

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>

#include <algorithm>
#include <iterator>
#include <map>
#include <type_traits>
#include <vector>

#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/dpg/errortools.hh>

namespace Dune {

namespace Imp {
  // TODO: implement std::hash or operator< for EntitySeed, so that I can
  // use unordered_map or map.
  // Alternatively, look into using Dune::PersistentContainer or look into
  // how dune-subgrid attaches data to grid cells.
  template<class Grid>
  class ProjectionMap
  {
    public:
      using Entity = typename Grid::LeafGridView::template Codim<0>::Entity;
      using EntitySeed = typename Entity::EntitySeed;
      using IdSet = typename Grid::GlobalIdSet;
      using IdType = typename IdSet::IdType;

      using mapped_type = Eigen::VectorXd;

      ProjectionMap(const IdSet& idSet)
        : projections{}
        , idSet(idSet)
        {}

      inline bool contains(const Entity& entity) const {
        return projections.find(idSet.id(entity)) != projections.cend();
      }

      const mapped_type& at(const Entity& entity) const {
        return projections.at(idSet.id(entity));
      }

      void insert(const Entity& entity, mapped_type&& localProjection) {
        projections.insert_or_assign(idSet.id(entity),
                                     std::move(localProjection));
      }

      void erase(const Entity& entity) {
        projections.erase(idSet.id(entity));
      }

    private:
      std::unordered_map<IdType, Eigen::VectorXd> projections;
      const IdSet& idSet;
  };

  template<class Grid, class GlobalBasis>
  std::vector<double>
  createProjectionVector(
    const ProjectionMap<Grid>& projection,
    GlobalBasis& globalBasis)
  {
    std::vector<double> projectionVector(globalBasis.size());
    auto localView = globalBasis.localView();
    for(const auto& e : elements(globalBasis.gridView()))
    {
      localView.bind(e);
      const auto& localCoefficients = projection.at(e);
      iterateOverLocalIndices(
        localView,
        [&](size_t i, auto gi)
        {
          projectionVector[gi[0]] = localCoefficients(i);
        },
        [](size_t i){},
        [](size_t i, auto gi, double wi)
        { /* will be copied from neighboring element */ }
      );
    }
    return projectionVector;
  }
} // End namespace Imp.

template<class Element, class LocalBasis, class LocalFunction>
std::tuple<Eigen::VectorXd, double>
projectLocalFunctionToEntity(
    const Element& element,
    const LocalBasis& localBasis,
    const LocalFunction& function)
{
  // Otherwise, we need to change the value type of the Eigen vector and matrix.
  static_assert(std::is_same_v<typename LocalBasis::Traits::RangeType,
                               FieldVector<double, 1>>);

  constexpr int dim = Element::mydimension;
  const auto geometry = element.geometry();

  // Quadrature order * 2 for square of polynomial FE function.
  // Take some extra degree in quadrature rule to estimate residual.
  const unsigned int quadratureOrder = 2 * localBasis.order() + 10;

  const auto& quad = Dune::QuadratureRules<double, dim>::rule(element.type(),
                                                              quadratureOrder);

  Eigen::VectorXd f(quad.size());
  {
    auto fIt = f.begin();
    for (const auto& quadPoint : quad) {
      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quadPoint.position();

      // The transformed quadrature weight
      const double integrationWeight
          = quadPoint.weight() * geometry.integrationElement(quadPos);

      *fIt = function(quadPos) * integrationWeight;
      ++fIt;
    }
  }

  using EigenRowMatrix
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  EigenRowMatrix A(quad.size(), localBasis.size());
  {
    std::vector<FieldVector<double,1>> shapeFunctionValues;
    shapeFunctionValues.reserve(localBasis.size());
    auto row = A.rowwise().begin();

    for (const auto& quadPoint : quad) {
      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quadPoint.position();

      // The transformed quadrature weight
      const double integrationWeight
          = quadPoint.weight() * geometry.integrationElement(quadPos);

      localBasis.evaluateFunction(quadPos, shapeFunctionValues);

      auto bx = shapeFunctionValues.cbegin();
      for(auto& entry : *row) {
        entry = (*bx)[0] * integrationWeight;
        ++bx;
      }
      ++row;
    }
  }

  // Project by solving least squares problem ||Ax - f||.
  // TODO: calculate projection and residual in one go.
  EigenRowMatrix Acopy(A);
  Eigen::ColPivHouseholderQR<Eigen::Ref<EigenRowMatrix>> decomposition(A);
  const auto solution = decomposition.solve(f);
  const double squaredError = (Acopy * solution - f).squaredNorm();
  return {solution, squaredError};
}

/** Project function to grid, refining until given accuracy is reached.
 *
 * This function projects a given function onto a finite element space on
 * the given grid. It automatically refines the grid via Dörfler marking
 * where necessary to attain a given accuracy.
 */
// TODO: Restrict domain of Function to coincide with global basis' domain.
template<class Grid, class GlobalBasis, class Function>
std::vector<double>
adaptivelyProjectFunctionToGrid(
    Grid& grid,
    GlobalBasis& globalBasis,
    Function& function,
    double doerflerRatio,
    double maxError)
{
  using GridView = typename GlobalBasis::GridView;
  static_assert(std::is_same_v<typename Grid::LeafGridView, GridView>);
  using Entity = typename GridView::template Codim<0>::Entity;
  using EntitySeed = typename Entity::EntitySeed;

  const GridView gridView = globalBasis.gridView();
  auto localView = globalBasis.localView();
  auto localF = localFunction(Functions::makeGridViewFunction(function,
                                                              gridView));

  const double maxSquaredError = maxError * maxError;

  Imp::ProjectionMap<Grid> projection{grid.globalIdSet()};
  std::vector<std::tuple<EntitySeed, double>> errorEstimates;
  errorEstimates.reserve(gridView.size(0));
  {
    for(const auto& e : elements(gridView))
    {
      localF.bind(e);
      localView.bind(e);
      const auto& localBasis = localView.tree().finiteElement().localBasis();

      auto [localProjection, squaredElementError]
        = projectLocalFunctionToEntity(e, localBasis, localF);
      projection.insert(e, std::move(localProjection));
      errorEstimates.emplace_back(e.seed(), squaredElementError);
    }
  }
  while(true) {
    std::stable_sort(errorEstimates.begin(), errorEstimates.end(),
        [](const auto& a, const auto& b) {
            return std::get<1>(a) > std::get<1>(b);
        });
    auto [squaredError, endDoerfler] = Doerfler::partition(doerflerRatio,
        errorEstimates.begin(), errorEstimates.end());
    if(squaredError < maxSquaredError) {
      return Imp::createProjectionVector(projection, globalBasis);
    } else {
      // TODO: Where do I get the grid from? Do I have to pass it explicitly?
      Doerfler::mark(grid, errorEstimates.begin(), endDoerfler);
      grid.preAdapt();
      grid.adapt();
      grid.postAdapt();
      globalBasis.update(grid.leafGridView());
      // Remove refined entities from projection and errorEstimates.
      for(auto it = errorEstimates.cbegin(); it != endDoerfler; ++it) {
        const auto e = grid.entity(std::get<0>(*it));
        projection.erase(e);
      }
      // Move endDoerfler..end() to front, and later remove refined elements.
      endDoerfler = std::move(endDoerfler, errorEstimates.end(),
                              errorEstimates.begin());
      // Remove entities that were refined to assure hanging node conditions.
      // TODO: If we do this before moving unrefined elements to the front,
      //       we can potentially save a few moves.
      for(auto it = errorEstimates.begin();
          // We cannot directly check it != endDoerfler, because removing
          // the last element would make us skip over endDoerfler, thus
          // causing us to overrun the end of errorEstimates.
          std::distance(it, endDoerfler) > 0;
          ++it)
      {
        const auto e = grid.entity(std::get<0>(*it));
        if(!e.isLeaf()) {
          projection.erase(e);
          std::iter_swap(it, endDoerfler);
          --endDoerfler;
        }
      }
      // Actually remove elements from errorEstimates.
      errorEstimates.erase(endDoerfler, errorEstimates.end());
      errorEstimates.reserve(gridView.size(0));

      for(const auto& e : elements(gridView))
      {
        if(projection.contains(e)) {
          continue;
        }
        localF.bind(e);
        localView.bind(e);
        const auto& localBasis = localView.tree().finiteElement().localBasis();

        auto [localProjection, squaredElementError]
          = projectLocalFunctionToEntity(e, localBasis, localF);
        projection.insert(e, std::move(localProjection));
        errorEstimates.emplace_back(e.seed(), squaredElementError);
      }
    }
  }
  // Could not reach requested accuracy in maximal number of iterations.
  return Imp::createProjectionVector(projection, globalBasis);
}

} // end namespace Dune

#endif
