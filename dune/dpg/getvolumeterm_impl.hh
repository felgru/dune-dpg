// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_GETVOLUMETERM_IMPL_HH
#define DUNE_DPG_GETVOLUMETERM_IMPL_HH

#include <dune/functions/common/functionconcepts.hh>
#include <dune/istl/bvector.hh>

#include "quadrature.hh"
#include "localevaluation.hh"

namespace Dune {
namespace detail {

// Compute the source term for a single element
template <LinearIntegrationType integrationType,
          class Space,
          bool = is_RefinedFiniteElement
                 <Space>::value>
struct GetVolumeTerm_Impl
{
  using LocalViewTest = typename Space::LocalView;

  template<class VectorType,
           class FactorType,
           class DirectionType>
  static void getVolumeTerm (const LocalViewTest& localViewTest,
                             VectorType& elementVector,
                             size_t spaceOffset,
                             const unsigned int quadratureOrder,
                             const FactorType& factor,
                             const DirectionType& beta);
};


template <LinearIntegrationType integrationType,
          class Space>
struct GetVolumeTerm_Impl<integrationType, Space, false>
{
  using LocalViewTest = typename Space::LocalView;

  template<class VectorType,
           class FactorType,
           class DirectionType>
  static void getVolumeTerm (const LocalViewTest& localViewTest,
                             VectorType& elementVector,
                             size_t spaceOffset,
                             const unsigned int quadratureOrder,
                             const FactorType& factor,
                             const DirectionType& beta)
  {
    static_assert(models<Functions::Concept::
         Function<double(const Dune::FieldVector<double, 2>&)>, FactorType>(),
         "The factor passed to getVolumeTerm does not model the "
         "Function concept.");

    using TestSpace = typename LocalViewTest::GlobalBasis;

    // Get the grid element from the local FE basis view
    typedef typename LocalViewTest::Element Element;
    const Element& element = localViewTest.element();

    const int dim = Element::dimension;
    auto geometry = element.geometry();

    // Get set of shape functions for this element
    const auto& localFiniteElementTest = localViewTest.tree().finiteElement();

    const unsigned int nDofs(localFiniteElementTest.localBasis().size());

    // TODO: This does not work with transport elements, as we do not know
    //       the transport direction.
    typename detail::ChooseQuadrature<TestSpace, TestSpace, Element>::type quad
      = detail::ChooseQuadrature<TestSpace, TestSpace, Element>
        ::Quadrature(element, quadratureOrder, nullptr);


    for ( size_t pt=0, qsize=quad.size(); pt < qsize; pt++ ) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();
      // Global position of the current quadrature point
      const FieldVector<double,dim>& globalQuadPos
          = geometry.global(quadPos);

      // The multiplicative factor in the integral transformation formula
      const double integrationWeight = geometry.integrationElement(quadPos)
                                       * quad[pt].weight()
                                       * factor(globalQuadPos);

     constexpr auto evaluationType = (integrationType ==
                                      LinearIntegrationType::valueFunction)
                                      ? EvaluationType::value : EvaluationType::grad;

     std::vector<FieldVector<double,1> > shapeFunctionValues =
     detail::LocalFunctionEvaluation<dim, evaluationType,
                                        DomainOfIntegration::interior>()
                      (localFiniteElementTest,
                       quadPos,
                       geometry,
                       beta);

      for (size_t i=0; i<nDofs; i++)
      {
        elementVector[i+spaceOffset] += shapeFunctionValues[i] * integrationWeight;
      }
    }

  }
};

template <LinearIntegrationType integrationType,
          class Space>
struct GetVolumeTerm_Impl<integrationType, Space, true>
{
  using LocalViewTest = typename Space::LocalView;

  template<class VectorType,
           class FactorType,
           class DirectionType>
  static void getVolumeTerm (const LocalViewTest& localViewTest,
                             VectorType& elementVector,
                             size_t spaceOffset,
                             const unsigned int quadratureOrder,
                             const FactorType& factor,
                             const DirectionType& beta)
  {
    static_assert(models<Functions::Concept::
         Function<double(const Dune::FieldVector<double, 2>&)>, FactorType>(),
         "The factor passed to getVolumeTerm does not model the "
         "Function concept.");

    using TestSpace = typename LocalViewTest::GlobalBasis;

    // Get the grid element from the local FE basis view
    typedef typename LocalViewTest::Element Element;
    const Element& element = localViewTest.element();

    const int dim = Element::dimension;
    auto geometry = element.geometry();

    // Get set of shape functions for this element
    const auto& localFiniteElementTest = localViewTest.tree().finiteElement();

    // TODO: This does not work with transport elements, as we do not know
    //       the transport direction.
    typename detail::ChooseQuadrature<TestSpace, TestSpace, Element>::type quad
      = detail::ChooseQuadrature<TestSpace, TestSpace, Element>
        ::Quadrature(element, quadratureOrder, nullptr);


    const auto& referenceGrid
      = localViewTest.tree().refinedReferenceElement();
    auto referenceGridView = referenceGrid.leafGridView();

    const unsigned int subElementStride =
        (is_DGRefinedFiniteElement<Space>::value) ?
          localFiniteElementTest.localBasis().size() : 0;

    unsigned int subElementOffset = 0;
    unsigned int subElementIndex = 0;
    for(const auto& subElement : elements(referenceGridView)) {
      auto subGeometryInReferenceElement = subElement.geometry();
      for ( size_t pt=0, qsize=quad.size(); pt < qsize; pt++ ) {

        // Position of the current quadrature point in the reference element
        const FieldVector<double,dim>& quadPos = quad[pt].position();
        const FieldVector<double,dim>& globalQuadPos
            = geometry.global(subGeometryInReferenceElement.global(quadPos));

        // The multiplicative factor in the integral transformation formula
        const double weightedfunctionValue
          = factor(globalQuadPos)
          * geometry.integrationElement(subGeometryInReferenceElement
                                                        .global(quadPos))
          * subGeometryInReferenceElement.integrationElement(quadPos)
          * quad[pt].weight();

        ////////////////////
        // Test Functions //
        ////////////////////
        constexpr auto evaluationType = (integrationType ==
                            LinearIntegrationType::valueFunction)
                            ? EvaluationType::value : EvaluationType::grad;

        std::vector<FieldVector<double,1> > shapeFunctionValues =
            detail::LocalRefinedFunctionEvaluation
                    <dim, evaluationType, DomainOfIntegration::interior,
                     is_ContinuouslyRefinedFiniteElement<TestSpace>::value>()
                          (localFiniteElementTest,
                           subElementIndex,
                           quadPos,
                           geometry,
                           subGeometryInReferenceElement,
                           beta);

        for (size_t i=0, rhsSize=shapeFunctionValues.size(); i<rhsSize; i++)
          elementVector[i+spaceOffset+subElementOffset] += shapeFunctionValues[i] * weightedfunctionValue;

      }
      if(is_DGRefinedFiniteElement<TestSpace>::value)
        subElementOffset += subElementStride;
      subElementIndex++;
    }
  }
};

}} // end namespace Dune::detail

#endif // DUNE_DPG_GETVOLUMETERM_IMPL_HH
