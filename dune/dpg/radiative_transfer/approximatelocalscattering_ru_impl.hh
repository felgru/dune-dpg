namespace detail {

template <class TestSpace,
          class SolutionSpace>
struct GetLocalApproximateScattering<TestSpace, SolutionSpace, true, false>
{
using TestLocalView = typename TestSpace::LocalView;
using SolutionLocalView = typename SolutionSpace::LocalView;
using SolutionLocalIndexSet = typename SolutionSpace::LocalIndexSet;

template <class VectorType,
          class Element,
          class KernelApproximation>
inline static void interiorImpl(
    const TestLocalView& testLocalView,
    const SolutionLocalView& solutionLocalView,
    VectorType& localScattering,
    const SolutionLocalIndexSet& solutionLocalIndexSet,
    size_t globalSolutionSpaceOffset,
    // unsigned int quadratureOrder,
    const Element& element,
    const KernelApproximation& kernelApproximation,
    const std::vector<BlockVector<FieldVector<double,1> >>& x,
    size_t si)
{
  const int dim = Element::dimension;
  auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& testLocalFiniteElement = testLocalView.tree().finiteElement();
  const auto& solutionLocalFiniteElement = solutionLocalView.tree().finiteElement();

  /* TODO:
   * - Adapt quadrature also to the kernel k
   */
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

      // The multiplicative factor in the integral transformation formula
      const double integrationWeight
        = geometry.integrationElement(subGeometryInReferenceElement
                                                      .global(quadPos))
        * subGeometryInReferenceElement.integrationElement(quadPos)
        * quad[pt].weight();

      // Evaluate all test shape function values at this quadrature point
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

      const size_t numS = x.size();
      Eigen::VectorXd uValues(numS);
      {
        std::vector<FieldVector<double,1>> shapeFunctionValues;
        solutionLocalFiniteElement.localBasis().
            evaluateFunction(quadPos, shapeFunctionValues);
        for( size_t scatteringAngle=0;
             scatteringAngle<numS; ++scatteringAngle) {
          double uValue = 0; // in direction of scatteringAngle
          // Evaluate all shape function values at this point
          std::vector<FieldVector<double,1>> shapeFunctionValues;
          solutionLocalFiniteElement.localBasis().evaluateFunction(
              subGeometryInReferenceElement.global(quadPos),
              shapeFunctionValues);
          for (size_t j=0; j<shapeFunctionValues.size(); j++)
          {
            /* This assumes that solutionLocalIndexSets and
             * globalSolutionSpaceOffset don't change for different
             * scattering angles.
             */
            auto row =
                solutionLocalIndexSet.index(j)[0]
              + globalSolutionSpaceOffset;
            uValue += x[scatteringAngle][row] * shapeFunctionValues[j];
          }
          uValues(scatteringAngle) = uValue;
        }
      }
      kernelApproximation.applyToVector(uValues);

      const double factor = uValues(si) * integrationWeight;
      for (size_t i=0, lb_size=testLocalFiniteElement.localBasis().size();
           i<lb_size; i++)
        localScattering[i+subElementOffset]
            += factor * testShapeFunctionValues[i];
    }
    if(is_DGRefinedFiniteElement<TestSpace>::value)
      subElementOffset += subElementStride;
    subElementIndex++;
  }
}

};

} // end namespace detail
