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
           class LocalCoefficients>
  static void getLocalVector(const LocalViewTest& localViewTest,
                             Vector& elementVector,
                             size_t spaceOffset,
                             const unsigned int quadratureOrder,
                             const LocalCoefficients& localCoefficients);
};


template <LinearIntegrationType integrationType,
          class Space>
struct GetLocalLinearTermVector<integrationType, Space, false>
{
  using LocalViewTest = typename Space::LocalView;

  template<class Vector,
           class LocalCoefficients>
  static void getLocalVector(const LocalViewTest& localViewTest,
                             Vector& elementVector,
                             size_t spaceOffset,
                             const unsigned int quadratureOrder,
                             const LocalCoefficients& localCoefficients)
  {
    using TestSpace = typename LocalViewTest::GlobalBasis;

    // Get the grid element from the local FE basis view
    typedef typename LocalViewTest::Element Element;
    const Element element = localViewTest.element();

    constexpr int dim = Element::mydimension;
    const auto geometry = element.geometry();

    // Get set of shape functions for this element
    const auto& localFiniteElementTest = localViewTest.tree().finiteElement();

    const unsigned int nDofs(localFiniteElementTest.size());

    typename detail::ChooseQuadrature<TestSpace, TestSpace, Element>::type quad
      = detail::ChooseQuadrature<TestSpace, TestSpace, Element>
        ::Quadrature(element, quadratureOrder);


    for (const auto& quadPoint : quad) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quadPoint.position();

      // The multiplicative factor in the integral transformation formula
      const double integrationWeight = geometry.integrationElement(quadPos)
                                       * quadPoint.weight()
                                       * localCoefficients.localFactor()(quadPos);

      using FunctionEvaluator
        = detail::LocalLinearTermFunctionEvaluation<dim, integrationType>;

      const std::vector<FieldVector<double,1> > shapeFunctionValues =
        FunctionEvaluator::evaluate
                      (localFiniteElementTest,
                       quadPos,
                       geometry,
                       localCoefficients);

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
           class LocalCoefficients>
  static void getLocalVector(LocalViewTest& localViewTest,
                             Vector& elementVector,
                             size_t spaceOffset,
                             const unsigned int quadratureOrder,
                             const LocalCoefficients& localCoefficients)
  {
    using TestSpace = typename LocalViewTest::GlobalBasis;

    // Get the grid element from the local FE basis view
    typedef typename LocalViewTest::Element Element;
    const Element element = localViewTest.element();

    constexpr int dim = Element::mydimension;
    const auto geometry = element.geometry();

    typename detail::ChooseQuadrature<TestSpace, TestSpace, Element>::type quad
      = detail::ChooseQuadrature<TestSpace, TestSpace, Element>
        ::Quadrature(element, quadratureOrder);

    const auto referenceGridView =
        localViewTest.tree().refinedReferenceElementGridView();

    unsigned int subElementOffset = 0;
    unsigned int subElementIndex = 0;
    localViewTest.resetSubElements();
    for(const auto& subElement : elements(referenceGridView)) {
      localViewTest.bindSubElement(subElement);

      // Get set of shape functions for this subElement
      const auto& localFiniteElementTest = localViewTest.tree().finiteElement();

      const auto subGeometryInReferenceElement = subElement.geometry();
      for (const auto& quadPoint : quad) {

        // Position of the current quadrature point in the reference element
        const FieldVector<double,dim>& quadPos = quadPoint.position();
        const FieldVector<double,dim> elementQuadPos
            = subGeometryInReferenceElement.global(quadPos);

        // The multiplicative factor in the integral transformation formula
        const double weightedfunctionValue
          = localCoefficients.localFactor()(elementQuadPos)
          * geometry.integrationElement(elementQuadPos)
          * subGeometryInReferenceElement.integrationElement(quadPos)
          * quadPoint.weight();

        ////////////////////
        // Test Functions //
        ////////////////////
        using FunctionEvaluator
          = detail::LocalRefinedLinearTermFunctionEvaluation
                <dim, integrationType>;

        const std::vector<FieldVector<double,1> > shapeFunctionValues =
            FunctionEvaluator::template evaluate
                    <is_ContinuouslyRefinedFiniteElement<TestSpace>::value>
                          (localFiniteElementTest,
                           subElementIndex,
                           quadPos,
                           geometry,
                           subGeometryInReferenceElement,
                           localCoefficients);

        auto entry = elementVector.begin() + spaceOffset + subElementOffset;
        for(const auto& shapeFunctionValue : shapeFunctionValues) {
          *entry += shapeFunctionValue * weightedfunctionValue;
          ++entry;
        }
      }
      if(is_DGRefinedFiniteElement<TestSpace>::value)
        subElementOffset += localFiniteElementTest.size();
      subElementIndex++;
    }
  }
};

}} // end namespace Dune::detail

#endif // DUNE_DPG_LOCALLINEARTERM_IMPL_HH
