// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_FUNCTIONS_SUBGRIDPROJECTION_HH
#define DUNE_DPG_FUNCTIONS_SUBGRIDPROJECTION_HH

#include <list>
#include <numeric>
#include <utility>
#include <vector>

#include <boost/hana.hpp>

#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>
#include <dune/dpg/functions/localindexsetiteration.hh>
#include <dune/dpg/functions/refinementinterpolation.hh>
#include <dune/dpg/integralterm.hh>
#include <dune/dpg/quadrature.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>

#if DUNE_DPG_USE_LEAST_SQUARES_INSTEAD_OF_CHOLESKY
#  include <dune/dpg/leastsquares.hh>
#else
#  include <dune/dpg/cholesky.hh>
#endif

namespace Dune {

namespace detail {
  template<class GlobalBasis>
  struct InterpolateOnCellLocalFunction {
    using FiniteElement = typename GlobalBasis::LocalView::Tree::FiniteElement;
    using Domain = FieldVector<double, GlobalBasis::GridView::dimension>;
    using Range = FieldVector<double, 1>;

    InterpolateOnCellLocalFunction(
        const FiniteElement& finiteElement,
        const std::vector<Range>& elementData)
      : finiteElement(finiteElement),
        elementData(elementData)
    {}

    void evaluate(const Domain& x, Range& y) const {
      y = Range(0);

      auto&& localBasis = finiteElement.localBasis();

      shapeFunctionValues.resize(localBasis.size());
      localBasis.evaluateFunction(x, shapeFunctionValues);

      for(size_t i = 0; i < localBasis.size(); i++) {
        y += elementData[i] * shapeFunctionValues[i];
      }
    }

    const FiniteElement& finiteElement;
    const std::vector<Range>& elementData;
    mutable std::vector<Range> shapeFunctionValues;
  };

  template<class HostGridGlobalBasis, class SubGridGlobalBasis>
  struct InterpolateOnSubCellLocalFunction {
    using FiniteElement
        = typename HostGridGlobalBasis::LocalView::Tree::FiniteElement;
    using Domain = FieldVector<double,
                               HostGridGlobalBasis::GridView::dimension>;
    using Range = FieldVector<double, 1>;
    using SubGeometryInReferenceElement
        = typename SubGridGlobalBasis::LocalView::Tree::RefinementGrid
              ::template Codim<0>::Entity::Geometry;

    InterpolateOnSubCellLocalFunction(
        const FiniteElement& finiteElement,
        const SubGeometryInReferenceElement& subGeometryInReferenceElement,
        const std::vector<Range>& elementData)
      : finiteElement(finiteElement),
        subGeometryInReferenceElement(subGeometryInReferenceElement),
        elementData(elementData)
    {}

    void evaluate(const Domain& x, Range& y) const {
      y = Range(0);

      auto&& localBasis = finiteElement.localBasis();

      shapeFunctionValues.resize(localBasis.size());
      localBasis.evaluateFunction(subGeometryInReferenceElement.global(x),
                                  shapeFunctionValues);

      for(size_t i = 0; i < localBasis.size(); i++) {
        y += elementData[i] * shapeFunctionValues[i];
      }
    }

    const FiniteElement& finiteElement;
    const SubGeometryInReferenceElement& subGeometryInReferenceElement;
    const std::vector<Range>& elementData;
    mutable std::vector<Range> shapeFunctionValues;
  };

  template<int dim, class HostGridElement, class SubGridElement>
  AffineGeometry<double, dim, dim>
  hostInSubGridCellGeometry(const HostGridElement& hostGridElement,
      const SubGridElement& subGridElement)
  {
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
    const auto referenceElement
        = Dune::referenceElement<double, dim>(hostGridElement.type());
#else
    const auto& referenceElement
        = ReferenceElements<double, dim>::general(hostGridElement.type());
#endif
    const auto hostGridCellGeometry = hostGridElement.geometry();
    const auto subGridCellGeometry = subGridElement.geometry();
    const size_t numVertices = referenceElement.size(dim);
    std::vector<FieldVector<double, dim>> vertices(numVertices);
    for(size_t i = 0; i < numVertices; i++) {
      vertices[i] = subGridCellGeometry.local(
                      hostGridCellGeometry.global(
                        referenceElement.position(i, dim)));
    }
    return AffineGeometry<double, dim, dim>(referenceElement, vertices);
  }

  template<class Element, class CellData, class SubGridLocalView,
           class HostGridLocalView,
           typename std::enable_if<!is_RefinedFiniteElement<typename
              SubGridLocalView::GlobalBasis>::value>::type* = nullptr>
  void computeProjectionRhs(const Element& e,
      const CellData& cellData,
      const SubGridLocalView& subGridLocalView,
      HostGridLocalView& hostGridLocalView,
      BlockVector<FieldVector<double,1>>& projectionRhs)
  {
    static_assert(!is_RefinedFiniteElement<typename
              HostGridLocalView::GlobalBasis>::value,
              "computeProjectionRhs only defined for unrefined HostGrid basis");
    projectionRhs = 0;
    const int dim = Element::mydimension;
    const auto& hostGrid = hostGridLocalView.globalBasis().gridView().grid();

    for(const auto& hostCellData : cellData) {
      const auto eHost = hostGrid.entity(hostCellData.first);
      const std::vector<FieldVector<double, 1>>& hostCellCoefficients
          = hostCellData.second;
      const auto hostCellEmbedding = hostInSubGridCellGeometry<dim>(eHost, e);

      const auto quadratureOrder
          = subGridLocalView.tree().finiteElement().localBasis().order()
          + hostGridLocalView.tree().finiteElement().localBasis().order();
      typename detail::ChooseQuadrature<
        typename SubGridLocalView::GlobalBasis,
        typename HostGridLocalView::GlobalBasis, Element>::type quad
          = detail::ChooseQuadrature<
              typename SubGridLocalView::GlobalBasis,
              typename HostGridLocalView::GlobalBasis, Element>
            ::Quadrature(e, quadratureOrder);

      for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {
        // Position of the current quadrature point in the reference element
        const FieldVector<double, dim>& quadPos = quad[pt].position();
        // Global position of the current quadrature point
        const FieldVector<double, dim>& subGridQuadPos
              = hostCellEmbedding.global(quadPos);

        // The multiplicative factor in the integral transformation formula
        const double integrationWeight
            = eHost.geometry().integrationElement(quadPos)
            * quad[pt].weight();

        const auto& subGridLocalFiniteElement
            = subGridLocalView.tree().finiteElement();
        std::vector<FieldVector<double, 1> > subGridValues;
        subGridLocalFiniteElement.localBasis().evaluateFunction(subGridQuadPos,
                                                                subGridValues);

        const auto& hostLocalFiniteElement
            = hostGridLocalView.tree().finiteElement();
        std::vector<FieldVector<double, 1> > hostValues;
        hostLocalFiniteElement.localBasis().evaluateFunction(quadPos,
                                                             hostValues);
        const FieldVector<double, 1> hostValue
            = std::inner_product(hostCellCoefficients.cbegin(),
                hostCellCoefficients.cend(), hostValues.cbegin(),
                FieldVector<double, 1>{0.});

        // compute projectionRhs[i] = <v_{sg,i}, \sum_j c_j v_{hg,j}>
        for (unsigned int i=0, imax=subGridValues.size(); i<imax; i++) {
          projectionRhs[i] += subGridValues[i] * hostValue * integrationWeight;
        }
      }
    }
  }

  struct PointInTriangleTest
  {
    template<class Geometry>
    PointInTriangleTest(const Geometry& triangle)
    {
      assert(triangle.type().isTriangle());
      using Point = FieldVector<double, 2>;
      const Point corner0 = triangle.corner(0);
      const Point corner1 = triangle.corner(1);
      const Point corner2 = triangle.corner(2);

      m[0][0] = corner0[0] - corner2[0];
      m[1][0] = corner0[1] - corner2[1];
      m[0][1] = corner1[0] - corner2[0];
      m[1][1] = corner1[1] - corner2[1];
      m.invert();

      const Point tmp{ -corner2[0], -corner2[1] };
      m.mv(tmp, r);
    }

    bool containsPoint(const FieldVector<double, 2>& point) const
    {
      // Convert point to baricentric cooradinates and test
      // for positivity.
      FieldVector<double, 2> bc(r);
      m.umv(point, bc);

      return bc[0] > 0
          && bc[1] > 0
          && 1-bc[0]-bc[1] > 0;
    }

  private:
    FieldMatrix<double, 2, 2> m;
    FieldVector<double, 2> r;
  };

  template<class Element, class CellData, class SubGridLocalView,
           class HostGridLocalView,
           typename std::enable_if<is_RefinedFiniteElement<typename
              SubGridLocalView::GlobalBasis>::value>::type* = nullptr>
  void computeProjectionRhs(const Element& e,
      const CellData& cellData,
      const SubGridLocalView& subGridLocalView,
      HostGridLocalView& hostGridLocalView,
      BlockVector<FieldVector<double,1>>& projectionRhs)
  {
    static_assert(!is_RefinedFiniteElement<typename
              HostGridLocalView::GlobalBasis>::value,
              "computeProjectionRhs only defined for unrefined HostGrid basis");
    using SubGridSpace = typename SubGridLocalView::GlobalBasis;

    projectionRhs = 0;
    const int dim = Element::mydimension;
    const auto& hostGrid = hostGridLocalView.globalBasis().gridView().grid();

    const auto referenceGridView =
        subGridLocalView.tree().refinedReferenceElement().leafGridView();

    const auto& subGridLocalFiniteElement
        = subGridLocalView.tree().finiteElement();
    const unsigned int subElementStride =
        (is_DGRefinedFiniteElement<SubGridSpace>::value) ?
          subGridLocalFiniteElement.localBasis().size() : 0;

    unsigned int subElementOffset = 0;
    unsigned int subElementIndex = 0;
    for(const auto& subElement : elements(referenceGridView)) {
      const auto subGeometryInReferenceElement = subElement.geometry();
      const PointInTriangleTest
          subElementTriangle(subGeometryInReferenceElement);

      for(const auto& hostCellData : cellData) {
        const auto eHost = hostGrid.entity(hostCellData.first);
        const auto eHostGeometry = eHost.geometry();
        const auto hostCellEmbedding = hostInSubGridCellGeometry<dim>(eHost, e);
        // Check if eHost lies in subElement.
        if(!subElementTriangle.containsPoint(hostCellEmbedding.center()))
          continue;

        const std::vector<FieldVector<double, 1>>& hostCellCoefficients
            = hostCellData.second;

        const auto quadratureOrder
            = subGridLocalFiniteElement.localBasis().order()
            + hostGridLocalView.tree().finiteElement().localBasis().order();
        typename detail::ChooseQuadrature<
          typename SubGridLocalView::GlobalBasis,
          typename HostGridLocalView::GlobalBasis, Element>::type quad
            = detail::ChooseQuadrature<
                typename SubGridLocalView::GlobalBasis,
                typename HostGridLocalView::GlobalBasis, Element>
              ::Quadrature(e, quadratureOrder);

        for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {
          // Position of the current quadrature point in the reference element
          const FieldVector<double, dim>& quadPos = quad[pt].position();
          // Global position of the current quadrature point
          const FieldVector<double, dim>& subGridQuadPos
                = subGeometryInReferenceElement.local(
                      hostCellEmbedding.global(quadPos));

          // The multiplicative factor in the integral transformation formula
          const double integrationWeight
              = eHostGeometry.integrationElement(quadPos)
              * quad[pt].weight();

          std::vector<FieldVector<double, 1> > subGridValues;
          boost::hana::eval_if(
              is_ContinuouslyRefinedFiniteElement<SubGridSpace>{},
              [&](auto _)
              {
                subGridLocalFiniteElement.localBasis()
                    .evaluateFunction(subElementIndex,
                                      subGridQuadPos,
                                      subGridValues);
              },
              [&](auto _)
              {
                subGridLocalFiniteElement.localBasis()
                    .evaluateFunction(subGridQuadPos, subGridValues);
              });

          const auto& hostLocalFiniteElement
              = hostGridLocalView.tree().finiteElement();
          std::vector<FieldVector<double, 1> > hostValues;
          hostLocalFiniteElement.localBasis().evaluateFunction(quadPos,
                                                               hostValues);
          const FieldVector<double, 1> hostValue
              = std::inner_product(hostCellCoefficients.cbegin(),
                  hostCellCoefficients.cend(), hostValues.cbegin(),
                  FieldVector<double, 1>{0.});

          // compute projectionRhs[i+subElementOffset]
          //            = <v_{sg,i}, \sum_j c_j v_{hg,j}>
          auto p = projectionRhs.begin() + subElementOffset;
          for (unsigned int i=0, imax=subGridValues.size(); i<imax; i++) {
            *p += subGridValues[i] * hostValue * integrationWeight;
            ++p;
          }
        }
      }
      if(is_DGRefinedFiniteElement<SubGridSpace>::value)
        subElementOffset += subElementStride;
      subElementIndex++;
    }
  }

  template<class SubGridGlobalBasis, class HostGridGlobalBasis>
  std::vector<FieldVector<double, 1>>
  projectCellDataToSubGrid(
      const typename SubGridGlobalBasis::GridView::template Codim<0>::Entity& e,
      const SubGridGlobalBasis& subGridGlobalBasis,
      const HostGridGlobalBasis& hostGridGlobalBasis,
      std::vector<std::pair<typename HostGridGlobalBasis::
                                LocalView::Element::EntitySeed,
                            std::vector<FieldVector<double, 1>>>>& cellData)
  {
    const auto eHost = subGridGlobalBasis.gridView().grid()
                                .template getHostEntity<0>(e);
    if(cellData.size() == 1) {
      // host and subGrid cell are the same → projection is identity
      std::vector<FieldVector<double, 1>> localData
          = std::exchange(cellData[0].second, {});
      cellData.clear();
      return boost::hana::eval_if(
        std::is_same<changeGridView_t<HostGridGlobalBasis,
                                      typename SubGridGlobalBasis::GridView>,
                     SubGridGlobalBasis>{},
        [&](auto _)
        {
          return localData;
        },
        [&](auto _)
        {
          // Different global bases on host and sub grid.
          // → Interpolate from hostGridGlobalBasis to subGridGlobalBasis.
          //   (we assume the subGridGlobalBasis to be superset of
          //   the hostGridGlobalBasis on the same level.)
          auto subGridGlobalBasisNode = subGridGlobalBasis.nodeFactory()
              .node(Dune::TypeTree::hybridTreePath());
          subGridGlobalBasisNode.bind(e);
          auto hostGridGlobalBasisNode = hostGridGlobalBasis.nodeFactory()
              .node(Dune::TypeTree::hybridTreePath());
          hostGridGlobalBasisNode.bind(eHost);
          auto&& subGridLocalFiniteElement
              = subGridGlobalBasisNode.finiteElement();
          auto&& hostGridLocalFiniteElement
              = hostGridGlobalBasisNode.finiteElement();
          return boost::hana::eval_if(
            is_RefinedFiniteElement<SubGridGlobalBasis>{},
            [&](auto _)
            {
              static_assert(is_DGRefinedFiniteElement<SubGridGlobalBasis>{},
                "Interpolation not implemented for continuously refined"
                " finite elements!");
              auto subGridLocalView = subGridGlobalBasis.localView();
              subGridLocalView.bind(e);
              std::vector<FieldVector<double, 1>>
                  interpolatedLocalData(subGridLocalView.size());

              const auto referenceGridView =
                  subGridLocalView.tree().refinedReferenceElement()
                                         .leafGridView();

              const unsigned int subElementStride =
                  (is_DGRefinedFiniteElement<SubGridGlobalBasis>::value) ?
                    subGridLocalFiniteElement.localBasis().size() : 0;

              std::vector<FieldVector<double, 1>> interpolatedSubGridLocalData;
              unsigned int subElementOffset = 0;
              unsigned int subElementIndex = 0;
              for(const auto& subElement : elements(referenceGridView)) {
                const auto subGeometryInReferenceElement
                    = subElement.geometry();
                auto hostGridFunction
                  = detail::InterpolateOnSubCellLocalFunction
                      <HostGridGlobalBasis, SubGridGlobalBasis>(
                          hostGridLocalFiniteElement,
                          subGeometryInReferenceElement,
                          localData);
                subGridLocalFiniteElement.localInterpolation()
                  .interpolate(hostGridFunction, interpolatedSubGridLocalData);

                std::copy(interpolatedSubGridLocalData.cbegin(),
                          interpolatedSubGridLocalData.cend(),
                          interpolatedLocalData.begin() + subElementOffset);

                if(is_DGRefinedFiniteElement<SubGridGlobalBasis>::value)
                  subElementOffset += subElementStride;
                subElementIndex++;
              }
              return interpolatedLocalData;
            },
            [&](auto _)
            {
              auto hostGridFunction
                = detail::InterpolateOnCellLocalFunction<HostGridGlobalBasis>(
                    hostGridLocalFiniteElement,
                    localData);
              std::vector<FieldVector<double, 1>> interpolatedLocalData;
              subGridLocalFiniteElement.localInterpolation()
                .interpolate(hostGridFunction, interpolatedLocalData);
              return interpolatedLocalData;
            }
          );
        }
      );
    } else if(cellData.size() > 1) {
      // project cellData to e
      auto subGridLocalView = subGridGlobalBasis.localView();
      subGridLocalView.bind(e);
      auto hostGridLocalView = hostGridGlobalBasis.localView();
      hostGridLocalView.bind(eHost);
      // create matrix and rhs of projection problem and solve it
      Matrix<FieldMatrix<double,1,1>>
          projectionMatrix(subGridLocalView.size(), subGridLocalView.size());
      {
        IntegralTerm<IntegrationType::valueValue,
             DomainOfIntegration::interior, double> integralTerm(1.0);
        integralTerm.getLocalMatrix(subGridLocalView, subGridLocalView,
            projectionMatrix, 0, 0);
      }
      BlockVector<FieldVector<double,1>>
          projectionRhs(subGridLocalView.size());
      computeProjectionRhs(e, cellData, subGridLocalView, hostGridLocalView,
          projectionRhs);
#if DUNE_DPG_USE_LEAST_SQUARES_INSTEAD_OF_CHOLESKY
      solveLeastSquares(projectionMatrix, projectionRhs);
#else
      {
        Cholesky<Matrix<FieldMatrix<double,1,1>>> cholesky(projectionMatrix);
        cholesky.apply(projectionRhs);
      }
#endif
      std::vector<FieldVector<double, 1>> projection(subGridLocalView.size());
      for(size_t i = 0; i < projection.size(); i++) {
        projection[i] = projectionRhs[i];
      }
      return projection;
    } else {
      DUNE_THROW(InvalidStateException,
                 "cellData is expected to be non-empty.");
    }
  }
} // end namespace detail

template<class SubGridGlobalBasis, class HostGridGlobalBasis, class Vector>
class SubGridProjectionData
{
  using SubGridEntitySeed
      = typename SubGridGlobalBasis::GridView::template Codim<0>::
                    Entity::EntitySeed;
  using HostGridEntitySeed
      = typename HostGridGlobalBasis::LocalView::Element::EntitySeed;

  using CellData = std::vector<std::pair<HostGridEntitySeed,
                               std::vector<FieldVector<double, 1>>>>;
  using GridData = std::list<std::tuple<
                                SubGridEntitySeed,
                                std::vector<FieldVector<double, 1>>,
                                CellData>>;

public:

  SubGridProjectionData(
      const SubGridGlobalBasis& subGridGlobalBasis,
      const HostGridGlobalBasis& hostGridGlobalBasis,
      const Vector& hostGridData)
  {
    assert(hostGridData.size() >= hostGridGlobalBasis.size());
    static_assert(std::is_same<typename HostGridGlobalBasis::GridView,
        typename HostGridGlobalBasis::GridView::Grid::LeafGridView>::value,
        "The HostGridGlobalBasis has to be defined on a LeafGridView!");
    auto localView = hostGridGlobalBasis.localView();
    auto localIndexSet = hostGridGlobalBasis.localIndexSet();

    const unsigned int maxHostGridLevel =
        hostGridGlobalBasis.gridView().grid().maxLevel();
    const auto& subGrid = subGridGlobalBasis.gridView().grid();
    for(const auto& e : elements(subGridGlobalBasis.gridView()))
    {
      CellData cellData;
      const auto eHost = subGrid.template getHostEntity<0>(e);
      if(eHost.isLeaf()) {
        localView.bind(eHost);
        localIndexSet.bind(localView);

        std::vector<FieldVector<double, 1>>
            hostGridLocalData(localIndexSet.size());
        iterateOverLocalIndexSet(
          localIndexSet,
          [&](size_t i, auto gi)
          {
            hostGridLocalData[i] = hostGridData[gi[0]];
          },
          [&](size_t i){ hostGridLocalData[i] = 0; },
          [&](size_t i, auto gi, double wi)
          {
            hostGridLocalData[i] += wi * hostGridData[gi[0]];
          }
        );
        cellData.reserve(1);
        // direct transfer of hostGridData
        // Will be later interpolated in the call to
        // projectCellDataToSubGrid if global bases on host and
        // sub grid differ.
        cellData.push_back(std::make_pair(eHost.seed(),
                                          std::move(hostGridLocalData)));
      } else { // e is not contained in the HostLeafGridView:
        for(const auto& child : descendantElements(eHost, maxHostGridLevel))
        {
          if(child.isLeaf())
          {
            localView.bind(child);
            localIndexSet.bind(localView);

            std::vector<FieldVector<double, 1>> hostGridLocalData(localIndexSet.size());
            iterateOverLocalIndexSet(
              localIndexSet,
              [&](size_t i, auto gi)
              {
                hostGridLocalData[i] = hostGridData[gi[0]];
              },
              [&](size_t i){ hostGridLocalData[i] = 0; },
              [&](size_t i, auto gi, double wi)
              {
                hostGridLocalData[i] += wi * hostGridData[gi[0]];
              }
            );
            cellData.push_back(std::make_pair(child.seed(),
                                              std::move(hostGridLocalData)));
          }
        }
      }

      std::vector<FieldVector<double, 1>> cellProjection
          = detail::projectCellDataToSubGrid(e, subGridGlobalBasis,
                                             hostGridGlobalBasis, cellData);
      gridData.push_back(std::make_tuple(e.seed(), cellProjection, cellData));
    }
  }

  /**
   * project/interpolate data saved with attachDataToSubGrid to refined grid
   *
   * \note This function assumes that the grid only includes one entity type.
   */
  std::vector<FieldVector<double, 1>>
  restoreDataToRefinedSubGrid(
      const SubGridGlobalBasis& subGridGlobalBasis)
  {
    using SubGridElement = typename SubGridGlobalBasis::LocalView::Element;

    auto subGridView = subGridGlobalBasis.gridView();
    auto& subGrid = subGridView.grid();
    auto localView = subGridGlobalBasis.localView();
    auto localIndexSet = subGridGlobalBasis.localIndexSet();
    auto node = subGridGlobalBasis.nodeFactory()
                  .node(Dune::TypeTree::hybridTreePath());

    HostGridGlobalBasis hostGridGlobalBasis(
        subGrid.getHostGrid().leafGridView());

    std::vector<FieldVector<double, 1>> data(subGridGlobalBasis.size());

    // Iterate over gridData. If cell is not leaf in subGrid project
    // or interpolate data to its children.
    for(auto currentData = gridData.begin(); currentData != gridData.end(); )
    {
      const auto e = subGrid.entity(std::get<0>(*currentData));
      localView.bind(e);
      localIndexSet.bind(localView);

      if(e.isLeaf()) {
        ++currentData;
      } else {
        // e has been refined
        if(std::get<2>(*currentData).size() > 1) {
          // project from hostGrid to children
          for (const auto& child : descendantElements(e, subGrid.maxLevel())) {
            CellData childData;
            const auto childLevel = child.level();
            const auto childInHostGrid
                = subGrid.template getHostEntity<0>(child);
            for(auto& hostLeafData : std::get<2>(*currentData)) {
              auto hostCell = subGrid.getHostGrid().entity(hostLeafData.first);
              for(auto level = hostCell.level(); level > childLevel; level--)
                hostCell = hostCell.father();
              if(hostCell == childInHostGrid) {
                childData.push_back(std::move(hostLeafData));
              }
            }
            std::vector<FieldVector<double, 1>> childProjection
              = detail::projectCellDataToSubGrid(child,
                                                 subGridGlobalBasis,
                                                 hostGridGlobalBasis,
                                                 childData);
            gridData.insert(currentData,
                std::make_tuple(child.seed(), childProjection, childData));
          }
          currentData = gridData.erase(currentData);
        } else { // if(std::get<2>(*currentData).size() == 0) {
          // refined beyond hostGrid → interpolate cell data to children
          using LocalData = std::vector<FieldVector<double, 1>>;
          const LocalData& localData = std::get<1>(*currentData);
          for (const auto& child : descendantElements(e, subGrid.maxLevel()))
          {
            node.bind(child);
            localView.bind(child);
            localIndexSet.bind(localView);

            // This assumes that e and child share the same finite element
            // and thus the same entity type.
            auto&& localFiniteElement = node.finiteElement();
            if(child.father() != e) {
              std::cerr << "e is not father of child!\n"
                << "e.level()=" << e.level()
                << "child.level()=" << child.level() << '\n';
              std::exit(1);
            }

            boost::hana::eval_if(
              is_RefinedFiniteElement<SubGridGlobalBasis>{},
              [&](auto _)
              {
                static_assert(is_DGRefinedFiniteElement<SubGridGlobalBasis>{},
                  "Interpolation not implemented for continuously refined"
                  " finite elements!");
                // With refinement level > 1 the embeddings get a lot more
                // complicated.
                static_assert(levelOfFE<SubGridGlobalBasis>::value <= 1,
                  "Interpolation only implemented for up to one level of"
                  " local refinement!");
                std::vector<FieldVector<double, 1>>
                    childLocalData(localView.size());

                const auto referenceGridView =
                    localView.tree().refinedReferenceElement()
                                           .leafGridView();

                const unsigned int subElementStride =
                    (is_DGRefinedFiniteElement<SubGridGlobalBasis>::value) ?
                      localFiniteElement.localBasis().size() : 0;

                using SubElement
                    = typename SubGridGlobalBasis::LocalView::Tree
                          ::RefinementGrid::template Codim<0>::Entity;
                using SubGeometryInReferenceElement
                    = typename SubElement::Geometry;

                constexpr int dim = 2;

                unsigned int sourceSubElementOffset = 0;
                SubElement sourceSubElement;
                for(const auto& sourceSubElement_
                    : elements(referenceGridView)) {
                  const SubGeometryInReferenceElement
                    sourceSubGeometryInReferenceElement
                      = sourceSubElement_.geometry();

                  const detail::PointInTriangleTest
                      subElementTriangle(sourceSubGeometryInReferenceElement);

                  const auto childEmbedding
                      = detail::hostInSubGridCellGeometry<dim>(child, e);
                  // Check if child lies in sourceSubElement.
                  if(subElementTriangle
                        .containsPoint(childEmbedding.center()))
                  {
                    sourceSubElement = sourceSubElement_;
                    break;
                  }

                  if(is_DGRefinedFiniteElement<SubGridGlobalBasis>::value)
                    sourceSubElementOffset += subElementStride;
                }
                assert(sourceSubElement != SubElement{});

                const SubGeometryInReferenceElement
                  sourceSubGeometryInReferenceElement
                    = sourceSubElement.geometry();

                unsigned int targetSubElementOffset = 0;
                for(const auto& targetSubElement
                    : elements(referenceGridView)) {
                  const auto targetSubGeometryInReferenceElement
                      = targetSubElement.geometry();

                  // Compute a Geometry that transformes from the
                  // target subElement to the source subElement.
                  const auto& geometryInFather = child.geometryInFather();
                  using SubGeometry = AffineGeometry<double, dim, dim>;
                  const SubGeometry subGeometry
                      ( child.type()
                      , sourceSubGeometryInReferenceElement.local(
                          geometryInFather.global(
                          targetSubGeometryInReferenceElement
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
                            .global(referenceElement<double, dim>
                              (child.type()).position(0,dim))))
#else
                            .global(ReferenceElements<double, dim>
                              ::general(child.type()).position(0,dim))))
#endif
                      , sourceSubGeometryInReferenceElement
                          .jacobianInverseTransposed({}).leftmultiply(
                            geometryInFather
                              .jacobianTransposed({}).leftmultiply(
                                targetSubGeometryInReferenceElement
                                  .jacobianTransposed({})))
                      );

                  const LocalData localDataSegment(
                      localData.cbegin() + sourceSubElementOffset,
                      localData.cbegin() + sourceSubElementOffset
                                         + localFiniteElement.size());
                  auto oldGridFunction
                    = detail::RestoreDataToRefinedGridFunction
                        <SubGridGlobalBasis,
                         SubGeometry,
                         const LocalData>(
                            localFiniteElement,
                            subGeometry,
                            localDataSegment);
                  std::vector<FieldVector<double, 1>> interpolatedData;
                  localFiniteElement.localInterpolation().interpolate(
                      oldGridFunction,
                      interpolatedData);

                  std::copy(interpolatedData.cbegin(),
                            interpolatedData.cend(),
                            childLocalData.begin() + targetSubElementOffset);

                  if(is_DGRefinedFiniteElement<SubGridGlobalBasis>::value)
                    targetSubElementOffset += subElementStride;
                }

                gridData.insert(currentData,
                    std::make_tuple(child.seed(), std::move(childLocalData),
                                    CellData{}));
              },
              [&](auto _)
              {
                auto oldGridFunction = detail::RestoreDataToRefinedGridFunction
                  <SubGridGlobalBasis,
                   typename SubGridElement::LocalGeometry,
                   LocalData>(
                      localFiniteElement,
                      child.geometryInFather(),
                      localData);
                std::vector<FieldVector<double, 1>> childLocalData;
                localFiniteElement.localInterpolation().interpolate(
                    oldGridFunction,
                    childLocalData);

                gridData.insert(currentData,
                    std::make_tuple(child.seed(), std::move(childLocalData),
                                    CellData{}));
              }
            );
          }
          currentData = gridData.erase(currentData);
        }
      }
    }

    // Iterate over gridData and transfer data to result vector.
    for(auto& currentData : gridData)
    {
      auto e = subGrid.entity(std::get<0>(currentData));
      localView.bind(e);
      localIndexSet.bind(localView);

      assert(e.isLeaf());
      // directly copy cell data
      auto& localData = std::get<1>(currentData);
      iterateOverLocalIndexSet(
        localIndexSet,
        [&](size_t i, auto gi)
        {
          data[gi[0]] = localData[i];
        },
        [](size_t i){},
        [](size_t i, auto gi, double wi)
        { /* will be copied from neighboring element */ }
      );
    }

    return data;
  }

private:
  GridData gridData;
};

/** create an object to project data on the host grid to a sub grid
 *
 * This function returns an object that stores for each element in
 * the \p subGridGlobalBasis the descendant elements in the host grid
 * together with the local data from \p hostGridData.
 *
 * The returned object can then later be used to get the projected data
 * after refinement of the sub grid by using its
 * SubGridProjectionData::restoreDataToRefinedSubGrid method.
 */
template<class SubGridGlobalBasis, class HostGridGlobalBasis, class Vector>
SubGridProjectionData<SubGridGlobalBasis, HostGridGlobalBasis, Vector>
attachDataToSubGrid(
    const SubGridGlobalBasis& subGridGlobalBasis,
    const HostGridGlobalBasis& hostGridGlobalBasis,
    const Vector& hostGridData)
{
  SubGridProjectionData<SubGridGlobalBasis, HostGridGlobalBasis, Vector>
    subGridData(subGridGlobalBasis, hostGridGlobalBasis, hostGridData);
  return subGridData;
}

} // end namespace Dune

#endif
