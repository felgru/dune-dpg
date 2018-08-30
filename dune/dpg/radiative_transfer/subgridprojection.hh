// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_FUNCTIONS_SUBGRIDPROJECTION_HH
#define DUNE_DPG_FUNCTIONS_SUBGRIDPROJECTION_HH

#include <algorithm>
#include <list>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/hana.hpp>

#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>
#include <dune/dpg/functions/gridviewfunctions.hh>
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
    static_assert(!is_RefinedFiniteElement<GlobalBasis>::value,
        "Interpolation does not work with refined GlobalBasis!");
    using FiniteElement = typename GlobalBasis::LocalView::Tree::FiniteElement;
    using Domain = FieldVector<double, GlobalBasis::GridView::dimension>;
    using Range = FieldVector<double, 1>;

    struct Traits
    {
       using DomainType = Domain;
       using RangeType  = Range;
    };

    InterpolateOnCellLocalFunction(
        const FiniteElement& finiteElement,
        const std::vector<Range>& elementData)
      : finiteElement(finiteElement),
        elementData(elementData)
    {}

    void evaluate(const Domain& x, Range& y) const {
      auto&& localBasis = finiteElement.localBasis();

      shapeFunctionValues.resize(localBasis.size());
      localBasis.evaluateFunction(x, shapeFunctionValues);

      y = std::inner_product(shapeFunctionValues.cbegin(),
                             shapeFunctionValues.cend(),
                             elementData.cbegin(),
                             Range(0));
    }

    const FiniteElement& finiteElement;
    const std::vector<Range>& elementData;
    mutable std::vector<Range> shapeFunctionValues;
  };

  template<class HostGridGlobalBasis, class SubGridGlobalBasis>
  struct InterpolateOnSubCellLocalFunction {
    static_assert(!is_RefinedFiniteElement<HostGridGlobalBasis>::value,
        "Interpolation does not work with refined HostGridGlobalBasis!");
    using FiniteElement
        = typename HostGridGlobalBasis::LocalView::Tree::FiniteElement;
    using Domain = FieldVector<double,
                               HostGridGlobalBasis::GridView::dimension>;
    using Range = FieldVector<double, 1>;
    using SubGeometryInReferenceElement
        = typename SubGridGlobalBasis::LocalView::Tree::RefinementGrid
              ::template Codim<0>::Entity::Geometry;

    struct Traits
    {
       using DomainType = Domain;
       using RangeType  = Range;
    };

    InterpolateOnSubCellLocalFunction(
        const FiniteElement& finiteElement,
        const SubGeometryInReferenceElement& subGeometryInReferenceElement,
        const std::vector<Range>& elementData)
      : finiteElement(finiteElement),
        subGeometryInReferenceElement(subGeometryInReferenceElement),
        elementData(elementData)
    {}

    void evaluate(const Domain& x, Range& y) const {
      auto&& localBasis = finiteElement.localBasis();

      shapeFunctionValues.resize(localBasis.size());
      localBasis.evaluateFunction(subGeometryInReferenceElement.global(x),
                                  shapeFunctionValues);

      y = std::inner_product(shapeFunctionValues.cbegin(),
                             shapeFunctionValues.cend(),
                             elementData.cbegin(),
                             Range(0));
    }

    const FiniteElement& finiteElement;
    const SubGeometryInReferenceElement& subGeometryInReferenceElement;
    const std::vector<Range>& elementData;
    mutable std::vector<Range> shapeFunctionValues;
  };

  template<class HostGridElement, class SubGridElement>
  AffineGeometry<typename HostGridElement::Geometry::ctype,
                 HostGridElement::mydimension,
                 HostGridElement::mydimension>
  hostInSubGridCellGeometry(const HostGridElement& hostGridElement,
      const SubGridElement& subGridElement)
  {
    static_assert(static_cast<int>(HostGridElement::mydimension)
                  == static_cast<int>(SubGridElement::mydimension),
        "HostGridElement and SubGridElement have different mydimension!");
    static_assert(std::is_same<typename HostGridElement::Geometry::ctype,
                               typename SubGridElement::Geometry::ctype>::value,
        "HostGridElement and SubGridElement have different ctype!");
    constexpr int dim = HostGridElement::mydimension;
    using ctype = typename HostGridElement::Geometry::ctype;
    const auto referenceElement
        = Dune::referenceElement<ctype, dim>(hostGridElement.type());
    const auto hostGridCellGeometry = hostGridElement.geometry();
    const auto subGridCellGeometry = subGridElement.geometry();
    const size_t numVertices = referenceElement.size(dim);
    std::vector<FieldVector<ctype, dim>> vertices(numVertices);
    for(size_t i = 0; i < numVertices; i++) {
      vertices[i] = subGridCellGeometry.local(
                      hostGridCellGeometry.global(
                        referenceElement.position(i, dim)));
    }
    return AffineGeometry<ctype, dim, dim>(referenceElement, vertices);
  }

  template<class Element, class CellData, class SubGridLocalView,
           class HostGridLocalView,
           typename std::enable_if<!is_RefinedFiniteElement<typename
              SubGridLocalView::GlobalBasis>::value>::type* = nullptr>
  void computeProjectionRhs(const Element& e,
      const CellData& cellData,
      SubGridLocalView& subGridLocalView,
      const HostGridLocalView& hostGridLocalView,
      BlockVector<FieldVector<double,1>>& projectionRhs)
  {
    static_assert(!is_RefinedFiniteElement<typename
              HostGridLocalView::GlobalBasis>::value,
              "computeProjectionRhs only defined for unrefined HostGrid basis");
    projectionRhs = 0;
    constexpr int dim = Element::mydimension;
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
      SubGridLocalView& subGridLocalView,
      const HostGridLocalView& hostGridLocalView,
      BlockVector<FieldVector<double,1>>& projectionRhs)
  {
    static_assert(!is_RefinedFiniteElement<typename
              HostGridLocalView::GlobalBasis>::value,
              "computeProjectionRhs only defined for unrefined HostGrid basis");
    using SubGridSpace = typename SubGridLocalView::GlobalBasis;

    projectionRhs = 0;
    constexpr int dim = Element::mydimension;
    const auto& hostGrid = hostGridLocalView.globalBasis().gridView().grid();

    const auto referenceGridView =
        subGridLocalView.tree().refinedReferenceElementGridView();

    unsigned int subElementOffset = 0;
    unsigned int subElementIndex = 0;
    subGridLocalView.resetSubElements();
    for(const auto& subElement : elements(referenceGridView)) {
      subGridLocalView.bindSubElement(subElement);
      const auto& subGridLocalFiniteElement
          = subGridLocalView.tree().finiteElement();

      const auto subGeometryInReferenceElement = subElement.geometry();
      const PointInTriangleTest
          subElementTriangle(subGeometryInReferenceElement);

      for(const auto& hostCellData : cellData) {
        const auto eHost = hostGrid.entity(hostCellData.first);
        const auto eHostGeometry = eHost.geometry();
        const auto hostCellEmbedding = hostInSubGridCellGeometry(eHost, e);
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
              [&](auto id)
              {
                id(subGridLocalFiniteElement).localBasis()
                    .evaluateFunction(subElementIndex,
                                      subGridQuadPos,
                                      subGridValues);
              },
              [&](auto id)
              {
                id(subGridLocalFiniteElement).localBasis()
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
        subElementOffset += subGridLocalFiniteElement.size();
      subElementIndex++;
    }
  }

  template<class SubGridGlobalBasis, class HostGridGlobalBasis,
           typename std::enable_if<
             is_RefinedFiniteElement<SubGridGlobalBasis>::value
           >::type* = nullptr>
  std::vector<FieldVector<double, 1>>
  interpolateDataOnSameCell(
      const typename SubGridGlobalBasis::GridView::template Codim<0>::Entity& e,
      const typename HostGridGlobalBasis::GridView::template Codim<0>::Entity& eHost,
      const SubGridGlobalBasis& subGridGlobalBasis,
      const HostGridGlobalBasis& hostGridGlobalBasis,
      const std::vector<FieldVector<double, 1>>& localData)
  {
    auto subGridLocalView = subGridGlobalBasis.localView();
    subGridLocalView.bind(e);
    auto hostGridLocalView = hostGridGlobalBasis.localView();
    hostGridLocalView.bind(eHost);
    auto&& hostGridLocalFiniteElement
        = hostGridLocalView.tree().finiteElement();
    static_assert(is_DGRefinedFiniteElement<SubGridGlobalBasis>::value,
      "Interpolation not implemented for continuously refined"
      " finite elements!");
    std::vector<FieldVector<double, 1>>
        interpolatedLocalData(subGridLocalView.size());
    auto interpolatedLocalDataIterator = interpolatedLocalData.begin();

    const auto referenceGridView =
        subGridLocalView.tree().refinedReferenceElementGridView();

    std::vector<FieldVector<double, 1>> interpolatedSubGridLocalData;
    subGridLocalView.resetSubElements();
    for(const auto& subElement : elements(referenceGridView)) {
      subGridLocalView.bindSubElement(subElement);
      auto&& subGridLocalFiniteElement
          = subGridLocalView.tree().finiteElement();

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

      interpolatedLocalDataIterator =
          std::copy(interpolatedSubGridLocalData.cbegin(),
                    interpolatedSubGridLocalData.cend(),
                    interpolatedLocalDataIterator);
    }
    return interpolatedLocalData;
  }

  template<class SubGridGlobalBasis, class HostGridGlobalBasis,
           typename std::enable_if<
             !is_RefinedFiniteElement<SubGridGlobalBasis>::value
           >::type* = nullptr>
  std::vector<FieldVector<double, 1>>
  interpolateDataOnSameCell(
      const typename SubGridGlobalBasis::GridView::template Codim<0>::Entity& e,
      const typename HostGridGlobalBasis::GridView::template Codim<0>::Entity& eHost,
      const SubGridGlobalBasis& subGridGlobalBasis,
      const HostGridGlobalBasis& hostGridGlobalBasis,
      const std::vector<FieldVector<double, 1>>& localData)
  {
    auto subGridLocalView = subGridGlobalBasis.localView();
    subGridLocalView.bind(e);
    auto hostGridLocalView = hostGridGlobalBasis.localView();
    hostGridLocalView.bind(eHost);
    auto&& subGridLocalFiniteElement
        = subGridLocalView.tree().finiteElement();
    auto&& hostGridLocalFiniteElement
        = hostGridLocalView.tree().finiteElement();
    auto hostGridFunction
      = detail::InterpolateOnCellLocalFunction<HostGridGlobalBasis>(
          hostGridLocalFiniteElement,
          localData);
    std::vector<FieldVector<double, 1>> interpolatedLocalData;
    subGridLocalFiniteElement.localInterpolation()
      .interpolate(hostGridFunction, interpolatedLocalData);
    return interpolatedLocalData;
  }

  template<class SubGridGlobalBasis, class HostGridGlobalBasis,
           typename std::enable_if<
             !std::is_same<changeGridView_t<HostGridGlobalBasis,
                                      typename SubGridGlobalBasis::GridView>,
                           SubGridGlobalBasis>::value>::type* = nullptr>
  std::vector<FieldVector<double, 1>>
  maybeInterpolateDataOnSameCell(
      const typename SubGridGlobalBasis::GridView::template Codim<0>::Entity& e,
      const typename HostGridGlobalBasis::GridView::template Codim<0>::Entity& eHost,
      const SubGridGlobalBasis& subGridGlobalBasis,
      const HostGridGlobalBasis& hostGridGlobalBasis,
      const std::vector<FieldVector<double, 1>>& localData)
  {
    // Different global bases on host and sub grid.
    // → Interpolate from hostGridGlobalBasis to subGridGlobalBasis.
    //   (we assume the subGridGlobalBasis to be a superset of
    //   the hostGridGlobalBasis on the same level.)
    return interpolateDataOnSameCell(e, eHost,
        subGridGlobalBasis, hostGridGlobalBasis, localData);
  }

  template<class SubGridGlobalBasis, class HostGridGlobalBasis,
           typename std::enable_if<
             std::is_same<changeGridView_t<HostGridGlobalBasis,
                                      typename SubGridGlobalBasis::GridView>,
                           SubGridGlobalBasis>::value>::type* = nullptr>
  std::vector<FieldVector<double, 1>>
  maybeInterpolateDataOnSameCell(
      const typename SubGridGlobalBasis::GridView::template Codim<0>::Entity& e,
      const typename HostGridGlobalBasis::GridView::template Codim<0>::Entity& eHost,
      const SubGridGlobalBasis& subGridGlobalBasis,
      const HostGridGlobalBasis& hostGridGlobalBasis,
      const std::vector<FieldVector<double, 1>>& localData)
  {
    return localData;
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
      return maybeInterpolateDataOnSameCell(e, eHost,
              subGridGlobalBasis, hostGridGlobalBasis, localData);
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
        auto oneFunc = Functions::makeConstantGridViewFunction(1.,
                                    subGridGlobalBasis.gridView());
        IntegralTerm<IntegrationType::valueValue,
                     DomainOfIntegration::interior,
                     detail::LocalCoefficients::OnlyFactor<decltype(oneFunc)>>
                                                        integralTerm(oneFunc);
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
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    auto localIndexSet = hostGridGlobalBasis.localIndexSet();
#endif

    const unsigned int maxHostGridLevel =
        hostGridGlobalBasis.gridView().grid().maxLevel();
    const auto subGridView = subGridGlobalBasis.gridView();
    const auto& subGrid = subGridView.grid();
    for(const auto& e : elements(subGridView))
    {
      CellData cellData;
      const auto eHost = subGrid.template getHostEntity<0>(e);
      if(eHost.isLeaf()) {
        localView.bind(eHost);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
        localIndexSet.bind(localView);
#endif

        std::vector<FieldVector<double, 1>>
            hostGridLocalData(localView.size());
        iterateOverLocalIndices(
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
          localView,
#else
          localIndexSet,
#endif
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
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
            localIndexSet.bind(localView);
#endif

            std::vector<FieldVector<double, 1>>
                hostGridLocalData(localView.size());
            iterateOverLocalIndices(
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
              localView,
#else
              localIndexSet,
#endif
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
      gridData.push_back(std::make_tuple(e.seed(),
                                         std::move(cellProjection),
                                         cellData));
    }
  }

  /**
   * project/interpolate data saved with attachDataToSubGrid to refined grid
   *
   * \note This function assumes that the grid only includes one entity type.
   */
  std::vector<FieldVector<double, 1>>
  restoreDataToRefinedSubGrid(
      const SubGridGlobalBasis& subGridGlobalBasis,
      // TODO: After refinement, the leaf grid of the hostGridGlobalBasis
      //       does not fit with the saved data anymore!
      const HostGridGlobalBasis& hostGridGlobalBasis)
  {
    transferDataToCellsOfRefinedSubGrid(subGridGlobalBasis,
                                        hostGridGlobalBasis);

    return createRefinedSubGridData(subGridGlobalBasis);
  }

private:

  void transferDataToCellsOfRefinedSubGrid(
      const SubGridGlobalBasis& subGridGlobalBasis,
      const HostGridGlobalBasis& hostGridGlobalBasis)
  {
    const auto subGridView = subGridGlobalBasis.gridView();
    const auto& subGrid = subGridView.grid();

    // Iterate over gridData. If cell is not leaf in subGrid project
    // or interpolate data to its children.
    for(auto currentData = gridData.begin(); currentData != gridData.end(); )
    {
      const auto e = subGrid.entity(std::get<0>(*currentData));
      if(e.isLeaf()) {
        ++currentData;
        continue;
      }

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
              std::make_tuple(child.seed(),
                std::move(childProjection), childData));
        }
        currentData = gridData.erase(currentData);
      } else /* if(std::get<2>(*currentData).size() == 0) */ {
        // refined beyond hostGrid → interpolate cell data to children
        auto sourceLocalView = subGridGlobalBasis.localView();
        sourceLocalView.bind(e);
        auto targetLocalView = subGridGlobalBasis.localView();
        using LocalData = std::vector<FieldVector<double, 1>>;
        const LocalData& sourceLocalData = std::get<1>(*currentData);
        for (const auto& child : descendantElements(e, subGrid.maxLevel()))
        {
          if(child.father() != e) {
            std::cerr << "e is not father of child!\n"
              << "e.level()=" << e.level()
              << "child.level()=" << child.level() << '\n';
            std::exit(1);
          }

          targetLocalView.bind(child);

          // This assumes that e and child share the same finite element
          // and thus the same entity type.
          boost::hana::eval_if(
            is_RefinedFiniteElement<SubGridGlobalBasis>{},
            [&](auto id)
            {
#if 0
              static_assert(is_DGRefinedFiniteElement<SubGridGlobalBasis>{},
                "Interpolation not implemented for continuously refined"
                " finite elements!");
#endif
              // With refinement level > 1 the embeddings get a lot more
              // complicated.
              static_assert(levelOfFE<SubGridGlobalBasis>::value <= 1,
                "Interpolation only implemented for up to one level of"
                " local refinement!");
              LocalData childLocalData(targetLocalView.size());

              // We assume that the referenceGridView of
              // source and target local view are the same.
              const auto referenceGridView =
                  id(sourceLocalView).tree().refinedReferenceElementGridView();

              using SubElement
                  = typename decltype(referenceGridView)
                        ::template Codim<0>::Entity;
              using SubGeometryInReferenceElement
                  = typename SubElement::Geometry;

              constexpr int dim = 2;

              auto sourceLocalDataBegin = sourceLocalData.cbegin();
              SubElement sourceSubElement;
              id(sourceLocalView).resetSubElements();
              for(const auto& sourceSubElement_ : elements(referenceGridView)) {
                id(sourceLocalView).bindSubElement(sourceSubElement_);
                auto&& sourceLocalFiniteElement
                    = sourceLocalView.tree().finiteElement();

                const SubGeometryInReferenceElement
                  sourceSubGeometryInReferenceElement
                    = sourceSubElement_.geometry();

                const detail::PointInTriangleTest
                    subElementTriangle(sourceSubGeometryInReferenceElement);

                const auto childEmbedding
                    = detail::hostInSubGridCellGeometry(child, e);
                // Check if child lies in sourceSubElement.
                if(subElementTriangle.containsPoint(childEmbedding.center()))
                {
                  sourceSubElement = sourceSubElement_;
                  break;
                }

                if(is_DGRefinedFiniteElement<SubGridGlobalBasis>::value)
                  sourceLocalDataBegin += sourceLocalFiniteElement.size();
              }
              assert(sourceSubElement != SubElement{});

              const SubGeometryInReferenceElement
                sourceSubGeometryInReferenceElement
                  = sourceSubElement.geometry();
              auto&& sourceLocalFiniteElement
                  = sourceLocalView.tree().finiteElement();

              auto childLocalDataIterator = childLocalData.begin();
              id(targetLocalView).resetSubElements();
              for(const auto& targetSubElement : elements(referenceGridView)) {
                id(targetLocalView).bindSubElement(targetSubElement);
                auto&& targetLocalFiniteElement
                    = targetLocalView.tree().finiteElement();

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
                          .global(referenceElement<double, dim>
                            (child.type()).position(0,dim))))
                    , sourceSubGeometryInReferenceElement
                        .jacobianInverseTransposed({}).leftmultiply(
                          geometryInFather
                            .jacobianTransposed({}).leftmultiply(
                              targetSubGeometryInReferenceElement
                                .jacobianTransposed({})))
                    );

                auto oldGridFunction
                  = detail::RestoreDataToRefinedGridFunction
                      <SubGridGlobalBasis,
                       SubGeometry,
                       decltype(sourceLocalDataBegin)>(
                          sourceLocalFiniteElement,
                          subGeometry,
                          sourceLocalDataBegin);
                std::vector<FieldVector<double, 1>> interpolatedData;
                targetLocalFiniteElement.localInterpolation().interpolate(
                    oldGridFunction,
                    interpolatedData);

                childLocalDataIterator =
                    std::copy(interpolatedData.cbegin(),
                              interpolatedData.cend(),
                              childLocalDataIterator);
              }

              gridData.insert(currentData,
                  std::make_tuple(child.seed(), std::move(childLocalData),
                                  CellData{}));
            },
            [&](auto id)
            {
              auto&& sourceLocalFiniteElement
                  = sourceLocalView.tree().finiteElement();
              auto&& targetLocalFiniteElement
                  = targetLocalView.tree().finiteElement();

              using SubGridElement
                  = typename SubGridGlobalBasis::LocalView::Element;
              auto oldGridFunction = detail::RestoreDataToRefinedGridFunction
                <SubGridGlobalBasis,
                 typename SubGridElement::LocalGeometry,
                 typename LocalData::const_iterator>(
                    sourceLocalFiniteElement,
                    id(child).geometryInFather(),
                    sourceLocalData.cbegin());
              std::vector<FieldVector<double, 1>> childLocalData;
              targetLocalFiniteElement.localInterpolation().interpolate(
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

  // This function assumes that gridData has already been transferred
  // to the refined subgrid by projection/interpolation.
  std::vector<FieldVector<double, 1>>
  createRefinedSubGridData(const SubGridGlobalBasis& subGridGlobalBasis) const
  {
    std::vector<FieldVector<double, 1>> data(subGridGlobalBasis.size());
    const auto& subGrid = subGridGlobalBasis.gridView().grid();
    auto localView = subGridGlobalBasis.localView();
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    auto localIndexSet = subGridGlobalBasis.localIndexSet();
#endif
    // Iterate over gridData and transfer data to result vector.
    for(auto& currentData : gridData)
    {
      auto e = subGrid.entity(std::get<0>(currentData));
      localView.bind(e);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
      localIndexSet.bind(localView);
#endif

      assert(e.isLeaf());
      // directly copy cell data
      auto& localData = std::get<1>(currentData);
      iterateOverLocalIndices(
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
        localView,
#else
        localIndexSet,
#endif
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
