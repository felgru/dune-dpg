// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_LOCALLINEARTERM_IMPL_HH
#define DUNE_DPG_LOCALLINEARTERM_IMPL_HH

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
struct GetLocalLinearTermVector
{
  using LocalViewTest = typename Space::LocalView;

  template<class Vector,
           class LocalFactor,
           class Direction>
  static void getLocalVector(const LocalViewTest& localViewTest,
                             Vector& elementVector,
                             size_t spaceOffset,
                             const unsigned int quadratureOrder,
                             const LocalFactor& localFactor,
                             const Direction& beta);
};


template <LinearIntegrationType integrationType,
          class Space>
struct GetLocalLinearTermVector<integrationType, Space, false>
{
  using LocalViewTest = typename Space::LocalView;

  template<class Vector,
           class LocalFactor,
           class Direction>
  static void getLocalVector(const LocalViewTest& localViewTest,
                             Vector& elementVector,
                             size_t spaceOffset,
                             const unsigned int quadratureOrder,
                             const LocalFactor& localFactor,
                             const Direction& beta)
  {
    static_assert(models<Functions::Concept::
         Function<double(const Dune::FieldVector<double, 2>&)>, LocalFactor>(),
         "The localFactor passed to getLocalVector does not model the "
         "Function concept.");

    using TestSpace = typename LocalViewTest::GlobalBasis;

    // Get the grid element from the local FE basis view
    typedef typename LocalViewTest::Element Element;
    const Element& element = localViewTest.element();

    constexpr int dim = Element::mydimension;
    const auto geometry = element.geometry();

    // Get set of shape functions for this element
    const auto& localFiniteElementTest = localViewTest.tree().finiteElement();

    const unsigned int nDofs(localFiniteElementTest.size());

    typename detail::ChooseQuadrature<TestSpace, TestSpace, Element>::type quad
      = detail::ChooseQuadrature<TestSpace, TestSpace, Element>
        ::Quadrature(element, quadratureOrder);


    for ( size_t pt=0, qsize=quad.size(); pt < qsize; pt++ ) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();

      // The multiplicative factor in the integral transformation formula
      const double integrationWeight = geometry.integrationElement(quadPos)
                                       * quad[pt].weight()
                                       * localFactor(quadPos);

     constexpr auto evaluationType = (integrationType ==
                                      LinearIntegrationType::valueFunction)
                                      ? EvaluationType::value : EvaluationType::grad;

     std::vector<FieldVector<double,1> > shapeFunctionValues =
     detail::LocalFunctionEvaluation<dim, evaluationType>()
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
struct GetLocalLinearTermVector<integrationType, Space, true>
{
  using LocalViewTest = typename Space::LocalView;

  template<class Vector,
           class LocalFactor,
           class Direction>
  static void getLocalVector(const LocalViewTest& localViewTest,
                             Vector& elementVector,
                             size_t spaceOffset,
                             const unsigned int quadratureOrder,
                             const LocalFactor& localFactor,
                             const Direction& beta)
  {
    static_assert(models<Functions::Concept::
         Function<double(const Dune::FieldVector<double, 2>&)>, LocalFactor>(),
         "The localFactor passed to getLocalVector does not model the "
         "Function concept.");

    using TestSpace = typename LocalViewTest::GlobalBasis;

    // Get the grid element from the local FE basis view
    typedef typename LocalViewTest::Element Element;
    const Element& element = localViewTest.element();

    constexpr int dim = Element::mydimension;
    const auto geometry = element.geometry();

    // Get set of shape functions for this element
    const auto& localFiniteElementTest = localViewTest.tree().finiteElement();

    typename detail::ChooseQuadrature<TestSpace, TestSpace, Element>::type quad
      = detail::ChooseQuadrature<TestSpace, TestSpace, Element>
        ::Quadrature(element, quadratureOrder);


    const auto referenceGridView =
        localViewTest.tree().refinedReferenceElementGridView();

    const unsigned int subElementStride =
        (is_DGRefinedFiniteElement<Space>::value) ?
          localFiniteElementTest.size() : 0;

    unsigned int subElementOffset = 0;
    unsigned int subElementIndex = 0;
    for(const auto& subElement : elements(referenceGridView)) {
      const auto subGeometryInReferenceElement = subElement.geometry();
      for ( size_t pt=0, qsize=quad.size(); pt < qsize; pt++ ) {

        // Position of the current quadrature point in the reference element
        const FieldVector<double,dim>& quadPos = quad[pt].position();
        const FieldVector<double,dim> elementQuadPos
            = subGeometryInReferenceElement.global(quadPos);

        // The multiplicative factor in the integral transformation formula
        const double weightedfunctionValue
          = localFactor(elementQuadPos)
          * geometry.integrationElement(elementQuadPos)
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
                    <dim, evaluationType,
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

#endif // DUNE_DPG_LOCALLINEARTERM_IMPL_HH
