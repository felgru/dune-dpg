namespace detail {

template <class TestSpace,
          class SolutionSpace>
struct GetLocalApproximateScattering<TestSpace, SolutionSpace, false, false>
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

  // Loop over all quadrature points
  for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPos);

    // Evaluate all test shape function values at this quadrature point
    std::vector<FieldVector<double,1>> testShapeFunctionValues;
    testLocalFiniteElement.localBasis()
        .evaluateFunction(quadPos, testShapeFunctionValues);

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
        solutionLocalFiniteElement.localBasis().
            evaluateFunction(quadPos, shapeFunctionValues);
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

    const double factor = uValues(si) * quad[pt].weight()
                          * integrationElement;
    for (size_t i=0, i_max=localScattering.size(); i<i_max; i++)
      localScattering[i] += factor * testShapeFunctionValues[i];
  }
}

};

} // end namespace detail
