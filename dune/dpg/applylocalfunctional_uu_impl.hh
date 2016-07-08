namespace detail {

template <class TestSpace,
          class SolutionSpace>
struct ApplyLocalFunctional<TestSpace, SolutionSpace, false, false>
{
using TestLocalView = typename TestSpace::LocalView;
using SolutionLocalView = typename SolutionSpace::LocalView;
using SolutionLocalIndexSet = typename SolutionSpace::LocalIndexSet;

template <class VectorType,
          class Element,
          class FunctionalVector>
inline static void interiorImpl(
    const TestLocalView& testLocalView,
    const SolutionLocalView& solutionLocalView,
    VectorType& elementVector,
    size_t spaceOffset,
    const SolutionLocalIndexSet& solutionLocalIndexSet,
    const Element& element,
    const FunctionalVector& functionalVector)
{
  const int dim = Element::dimension;
  auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& testLocalFiniteElement = testLocalView.tree().finiteElement();
  const auto& solutionLocalFiniteElement = solutionLocalView.tree().finiteElement();

  BlockVector<FieldVector<double,1>>
      localFunctionalVector(solutionLocalView.size());
  for (size_t j=0, jMax=localFunctionalVector.size(); j<jMax; j++)
  {
    auto row = solutionLocalIndexSet.index(j)[0];
    localFunctionalVector[j] = functionalVector[row];
  }

  const unsigned int quadratureOrder
      = solutionLocalFiniteElement.localBasis().order()
        + testLocalFiniteElement.localBasis().order();

  typename detail::ChooseQuadrature<TestSpace, SolutionSpace, Element>::type quad
    = detail::ChooseQuadrature<TestSpace, SolutionSpace, Element>
      ::Quadrature(element, quadratureOrder, nullptr);

  // Loop over all quadrature points
  for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPos);

    // Evaluate all shape function values at this quadrature point
    std::vector<FieldVector<double,1>> testShapeFunctionValues;
    testLocalFiniteElement.localBasis()
        .evaluateFunction(quadPos, testShapeFunctionValues);
    std::vector<FieldVector<double,1>> shapeFunctionValues;
    solutionLocalFiniteElement.localBasis().
        evaluateFunction(quadPos, shapeFunctionValues);

    double functionalValue = 0;
    for (size_t j=0, j_max=shapeFunctionValues.size(); j<j_max; j++)
      functionalValue += localFunctionalVector[j]
                       * shapeFunctionValues[j];
    functionalValue *= quad[pt].weight() * integrationElement;
    for (size_t i=0, i_max=elementVector.size(); i<i_max; i++) {
      elementVector[i+spaceOffset] += functionalValue * testShapeFunctionValues[i];
    }
  }
}

};

} // end namespace detail
