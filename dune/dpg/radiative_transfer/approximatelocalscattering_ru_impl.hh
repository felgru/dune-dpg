namespace detail {

template <class TestSpace,
          class SolutionSpace>
struct ApplyLocalFunctional<TestSpace, SolutionSpace, true, false>
{
using TestLocalView = typename TestSpace::LocalView;
using SolutionLocalView = typename SolutionSpace::LocalView;
using SolutionLocalIndexSet = typename SolutionSpace::LocalIndexSet;

template <class VectorType,
          class Element,
          class ScatteringFunctional>
inline static void interiorImpl(
    const TestLocalView& testLocalView,
    const SolutionLocalView& solutionLocalView,
    VectorType& localScattering,
    const SolutionLocalIndexSet& solutionLocalIndexSet,
    size_t globalSolutionSpaceOffset,
    const Element& element,
    const ScatteringFunctional& scatteringFunctional)
{
  const int dim = Element::dimension;
  auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& testLocalFiniteElement = testLocalView.tree().finiteElement();
  const auto& solutionLocalFiniteElement = solutionLocalView.tree().finiteElement();

  BlockVector<FieldVector<double,1>>
      localScatteringFunctional(solutionLocalView.size());
  for (size_t j=0, jMax=localScatteringFunctional.size(); j<jMax; j++)
  {
    auto row = solutionLocalIndexSet.index(j)[0];
    localScatteringFunctional[j] = scatteringFunctional[row];
  }

  const unsigned int quadratureOrder
      = solutionLocalFiniteElement.localBasis().order()
        + testLocalFiniteElement.localBasis().order();

  typename detail::ChooseQuadrature<TestSpace, SolutionSpace, Element>::type quad
    = detail::ChooseQuadrature<TestSpace, SolutionSpace, Element>
      ::Quadrature(element, quadratureOrder, nullptr);

  const auto& referenceGrid
    = testLocalView.tree().refinedReferenceElement();
  auto referenceGridView = referenceGrid.leafGridView();

  assert(element.type().isTriangle() || element.type().isQuadrilateral());
  const size_t subElementStride =
    (element.type().isTriangle())
    ? testLocalView.globalBasis().nodeFactory().dofsPerSubTriangle
    : testLocalView.globalBasis().nodeFactory().dofsPerSubQuad;

  unsigned int subElementOffset = 0;
  unsigned int subElementIndex = 0;
  for(const auto& subElement : elements(referenceGridView)) {
    auto subGeometryInReferenceElement = subElement.geometry();
    for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();
      const FieldVector<double,dim>& quadPosInReferenceElement =
          subGeometryInReferenceElement.global(quadPos);
      // The multiplicative factor in the integral transformation formula
      const double integrationWeight
        = geometry.integrationElement(quadPosInReferenceElement)
        * subGeometryInReferenceElement.integrationElement(quadPos)
        * quad[pt].weight();

      // Evaluate all shape function values at this quadrature point
      std::vector<FieldVector<double,1>> testShapeFunctionValues =
          detail::LocalRefinedFunctionEvaluation
                  <dim, EvaluationType::value, DomainOfIntegration::interior,
                   is_ContinuouslyRefinedFiniteElement<TestSpace>::value>()
                        (testLocalFiniteElement,
                         subElementIndex,
                         quadPos,
                         geometry,
                         subGeometryInReferenceElement,
                         {});
      std::vector<FieldVector<double,1>> shapeFunctionValues;
      solutionLocalFiniteElement.localBasis().
          evaluateFunction(quadPosInReferenceElement, shapeFunctionValues);

      double functionalValue = 0;
      for (size_t j=0, j_max=shapeFunctionValues.size(); j<j_max; j++)
        functionalValue += localScatteringFunctional[j]
                         * shapeFunctionValues[j];
      functionalValue *= integrationWeight;
      for (size_t i=0, i_max=testLocalFiniteElement.localBasis().size();
           i<i_max; i++) {
        localScattering[i+subElementOffset]
            += functionalValue * testShapeFunctionValues[i];
      }
    }
    if(is_DGRefinedFiniteElement<TestSpace>::value)
      subElementOffset += subElementStride;
    subElementIndex++;
  }
}

};

} // end namespace detail
