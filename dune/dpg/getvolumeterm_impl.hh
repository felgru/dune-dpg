// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_GETVOLUMETERM_IMPL_HH
#define DUNE_DPG_GETVOLUMETERM_IMPL_HH

#include <boost/fusion/sequence/intrinsic/at.hpp>

#include <dune/functions/common/functionconcepts.hh>
#include <dune/istl/bvector.hh>

#include "quadrature.hh"
#include "localevaluation.hh"

namespace Dune {
namespace detail {

// Compute the source term for a single element
template <class LocalViewTest, class VolumeTerm,
          bool = is_RefinedFiniteElement
                 <typename LocalViewTest::GlobalBasis>::value>
struct GetVolumeTerm_Impl
{
  void operator() (const LocalViewTest& localViewTest,
                   BlockVector<FieldVector<double,1> >& localRhs,
                   const VolumeTerm& volumeTerm);
};

template <class LocalViewTest, class VolumeTerm>
struct GetVolumeTerm_Impl<LocalViewTest, VolumeTerm, false>
{
  void operator() (const LocalViewTest& localViewTest,
                   BlockVector<FieldVector<double,1> >& localRhs,
                   const VolumeTerm& volumeTerm)
  {
    static_assert(models<Functions::Concept::
         Function<double(const Dune::FieldVector<double, 2>&)>, VolumeTerm>(),
         "The volumeTerm passed to getVolumeTerm does not model the "
         "Function concept.");

    using TestSpace = typename LocalViewTest::GlobalBasis;

    // Get the grid element from the local FE basis view
    typedef typename LocalViewTest::Element Element;
    const Element& element = localViewTest.element();

    const int dim = Element::dimension;
    auto geometry = element.geometry();

    // Get set of shape functions for this element
    const auto& localFiniteElementTest = localViewTest.tree().finiteElement();

    // Set all entries to zero
    localRhs.resize(localFiniteElementTest.localBasis().size());
    localRhs = 0;

    /* TODO: Quadrature order is only good enough for a constant volumeTerm. */
    const unsigned int quadratureOrder
        = localFiniteElementTest.localBasis().order();

    // TODO: This does not work with transport elements, as we do not know
    //       the transport direction.
    typename detail::ChooseQuadrature<TestSpace, TestSpace, Element>::type quad
      = detail::ChooseQuadrature<TestSpace, TestSpace, Element>
        ::Quadrature(element, quadratureOrder, nullptr);


    for ( size_t pt=0, qsize=quad.size(); pt < qsize; pt++ ) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();
      const FieldVector<double,dim>& globalQuadPos
          = geometry.global(quadPos);

      // The multiplicative factor in the integral transformation formula
      const double integrationElement
        = element.geometry().integrationElement(quadPos);

      const double weightedfunctionValue
        = volumeTerm(globalQuadPos) * quad[pt].weight() * integrationElement;

      std::vector<FieldVector<double,1> > shapeFunctionValues;
      localFiniteElementTest.localBasis().evaluateFunction(quadPos,
                                                           shapeFunctionValues);

      for (size_t i=0, rhsSize=localRhs.size(); i<rhsSize; i++)
        localRhs[i] += shapeFunctionValues[i] * weightedfunctionValue;

    }

  }
};

template <class LocalViewTest, class VolumeTerm>
struct GetVolumeTerm_Impl<LocalViewTest, VolumeTerm, true>
{
  void operator() (const LocalViewTest& localViewTest,
                   BlockVector<FieldVector<double,1> >& localRhs,
                   const VolumeTerm& volumeTerm)
  {
    static_assert(models<Functions::Concept::
         Function<double(const Dune::FieldVector<double, 2>&)>, VolumeTerm>(),
         "The volumeTerm passed to getVolumeTerm does not model the "
         "Function concept.");

    using TestSpace = typename LocalViewTest::GlobalBasis;

    // Get the grid element from the local FE basis view
    typedef typename LocalViewTest::Element Element;
    const Element& element = localViewTest.element();

    const int dim = Element::dimension;
    auto geometry = element.geometry();

    // Get set of shape functions for this element
    const auto& localFiniteElementTest = localViewTest.tree().finiteElement();

    // Set all entries to zero
    localRhs.resize(localFiniteElementTest.localBasis().size());
    localRhs = 0;

    /* TODO: Quadrature order is only good enough for a constant volumeTerm. */
    const unsigned int quadratureOrder
        = localFiniteElementTest.localBasis().order();

    // TODO: This does not work with transport elements, as we do not know
    //       the transport direction.
    typename detail::ChooseQuadrature<TestSpace, TestSpace, Element>::type quad
      = detail::ChooseQuadrature<TestSpace, TestSpace, Element>
        ::Quadrature(element, quadratureOrder, nullptr);


    const auto& referenceGrid
      = localViewTest.tree().refinedReferenceElement();
    auto referenceGridView = referenceGrid.leafGridView();

    assert(element.type().isTriangle() || element.type().isQuadrilateral());
    const size_t subElementStride =
      (element.type().isTriangle())
      ? localViewTest.globalBasis().nodeFactory().dofsPerSubTriangle
      : localViewTest.globalBasis().nodeFactory().dofsPerSubQuad;

    unsigned int subElementOffset = 0;
    unsigned int subElementIndex = 0;
    for(const auto& subElement : elements(referenceGridView)) {
      auto subGeometryInReferenceElement = subElement.geometry();
      for ( size_t pt=0, qsize=quad.size(); pt < qsize; pt++ ) {

        // Position of the current quadrature point in the reference element
        const FieldVector<double,dim>& quadPos = quad[pt].position();
        const FieldVector<double,dim>& globalQuadPos
            = geometry.global(subGeometryInReferenceElement.global(quadPos));

        // The transposed inverse Jacobian of the map from the reference element to the element
        const auto& jacobianSub
            = subGeometryInReferenceElement.jacobianInverseTransposed(quadPos);
        const auto& jacobian = geometry.jacobianInverseTransposed
                               (subGeometryInReferenceElement.global(quadPos));

        // The multiplicative factor in the integral transformation formula
        const double weightedfunctionValue
          = volumeTerm(globalQuadPos)
          * geometry.integrationElement(subGeometryInReferenceElement
                                                        .global(quadPos))
          * subGeometryInReferenceElement.integrationElement(quadPos)
          * quad[pt].weight();

        ////////////////////
        // Test Functions //
        ////////////////////
        std::vector<FieldVector<double,1> > shapeFunctionValues =
            detail::LocalRefinedFunctionEvaluation
                    <dim, EvaluationType::value, DomainOfIntegration::interior,
                     is_ContinuouslyRefinedFiniteElement<TestSpace>::value>()
                          (localFiniteElementTest,
                           subElementIndex,
                           quadPos,
                           geometry,
                           subGeometryInReferenceElement,
                           {});

        for (size_t i=0, rhsSize=localRhs.size(); i<rhsSize; i++)
          localRhs[i] += shapeFunctionValues[i] * weightedfunctionValue;

      }
      if(is_DGRefinedFiniteElement<TestSpace>::value)
        subElementOffset += subElementStride;
      subElementIndex++;
    }
  }
};

template <class LocalViewTest, class VolumeTerm>
inline void getVolumeTerm(const LocalViewTest& localViewTest,
                          BlockVector<FieldVector<double,1> >& localRhs,
                          const VolumeTerm& volumeTerm)
{
  GetVolumeTerm_Impl<LocalViewTest, VolumeTerm>()
          (localViewTest, localRhs, volumeTerm);
}

struct getVolumeTermHelper
{
  template<class Seq>
  void operator()(const Seq& seq) const
  {
    using namespace boost::fusion;

    // TODO: this can probably be done more elegantly by sequence fusion.
    getVolumeTerm(at_c<0>(seq), at_c<1>(seq), at_c<2>(seq));
  }
};

}} // end namespace Dune::detail

#endif // DUNE_DPG_GETVOLUMETERM_IMPL_HH
