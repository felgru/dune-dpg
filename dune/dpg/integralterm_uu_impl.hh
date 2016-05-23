namespace detail {

template <IntegrationType type,
          class LhsSpace,
          class RhsSpace>
struct GetLocalMatrix<type, LhsSpace, RhsSpace, false, false>
{
using LhsLocalView = typename LhsSpace::LocalView;
using RhsLocalView = typename RhsSpace::LocalView;

template <class MatrixType,
          class Element,
          class FactorType,
          class DirectionType>
inline static void interiorImpl(const LhsLocalView& lhsLocalView,
                                const RhsLocalView& rhsLocalView,
                                MatrixType& elementMatrix,
                                size_t lhsSpaceOffset,
                                size_t rhsSpaceOffset,
                                unsigned int quadratureOrder,
                                const Element& element,
                                const FactorType& factor,
                                const DirectionType& lhsBeta,
                                const DirectionType& rhsBeta)
{
  const int dim = Element::dimension;
  auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& lhsLocalFiniteElement = lhsLocalView.tree().finiteElement();
  const auto& rhsLocalFiniteElement = rhsLocalView.tree().finiteElement();

  const unsigned int nLhs(lhsLocalFiniteElement.localBasis().size());
  const unsigned int nRhs(rhsLocalFiniteElement.localBasis().size());

  typename detail::ChooseQuadrature<LhsSpace, RhsSpace, Element>::type quad
    = detail::ChooseQuadrature<LhsSpace, RhsSpace, Element>
      ::Quadrature(element, quadratureOrder, lhsBeta);

  for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();
    // Global position of the current quadrature point
    const FieldVector<double,dim>& globalQuadPos
          = geometry.global(quadPos);

    // The multiplicative factor in the integral transformation formula
    const double integrationWeight = geometry.integrationElement(quadPos)
                                   * quad[pt].weight()
                                   * detail::evaluateFactor(factor, globalQuadPos);


    //////////////////////////////
    // Left hand side Functions //
    //////////////////////////////
    constexpr auto lhsType = (type == IntegrationType::valueValue ||
                              type == IntegrationType::valueGrad)
                             ? EvaluationType::value : EvaluationType::grad;

    std::vector<FieldVector<double,1> > lhsValues =
        detail::LocalFunctionEvaluation<dim, lhsType,
                                        DomainOfIntegration::interior>()
                      (lhsLocalFiniteElement,
                       quadPos,
                       geometry,
                       lhsBeta);

    ///////////////////////////////
    // Right hand side Functions //
    ///////////////////////////////
    constexpr auto rhsType = (type == IntegrationType::valueValue ||
                              type == IntegrationType::gradValue)
                             ? EvaluationType::value : EvaluationType::grad;

    std::vector<FieldVector<double,1> > rhsValues =
        detail::LocalFunctionEvaluation<dim, rhsType,
                                        DomainOfIntegration::interior>()
                      (rhsLocalFiniteElement,
                       quadPos,
                       geometry,
                       rhsBeta);

    // Compute the actual matrix entries
    for (unsigned int i=0; i<nLhs; i++)
    {
      for (unsigned int j=0; j<nRhs; j++)
      {
        elementMatrix[i+lhsSpaceOffset][j+rhsSpaceOffset]
                += (lhsValues[i] * rhsValues[j]) * integrationWeight;
      }
    }
  }
}


template <class MatrixType,
          class Element,
          class FactorType,
          class DirectionType>
inline static void
faceImpl(const LhsLocalView& lhsLocalView,
         const RhsLocalView& rhsLocalView,
         MatrixType& elementMatrix,
         size_t lhsSpaceOffset,
         size_t rhsSpaceOffset,
         unsigned int quadratureOrder,
         const Element& element,
         const FactorType& factor,
         const DirectionType& lhsBeta,
         const DirectionType& rhsBeta)
{
  const int dim = Element::dimension;

  // Get set of shape functions for this element
  const auto& lhsLocalFiniteElement = lhsLocalView.tree().finiteElement();
  const auto& rhsLocalFiniteElement = rhsLocalView.tree().finiteElement();

  const unsigned int nLhs(lhsLocalFiniteElement.localBasis().size());
  const unsigned int nRhs(rhsLocalFiniteElement.localBasis().size());

  const auto& gridView = lhsLocalView.globalBasis().gridView();

  unsigned int nInflowIntersections = 0;
  unsigned int nOutflowIntersections = 0;
  for (auto&& intersection : intersections(gridView, element))
  {
    double prod = lhsBeta * intersection.centerUnitOuterNormal();
    if(prod > 0)
      ++nOutflowIntersections;
    else if (prod < 0)
      ++nInflowIntersections;
  }

  FieldVector<double,dim> referenceBeta;
  {
    const auto& jacobianInverse = element.geometry().jacobianInverseTransposed({0., 0.});
    jacobianInverse.mtv(lhsBeta, referenceBeta);
  }

  for (auto&& intersection : intersections(gridView, element))
  {
    if(type == IntegrationType::travelDistanceWeighted &&
       lhsBeta * intersection.centerUnitOuterNormal() >= 0) {
      /* Only integrate over inflow boundaries. */
      continue;
    }

    using intersectionType
      = typename std::decay<decltype(intersection)>::type;
    // TODO: Do we really want to have a transport quadrature rule
    //       on the faces, if one of the FE spaces is a transport space?
    QuadratureRule<double, 1> quadFace
      = detail::ChooseQuadrature<LhsSpace, RhsSpace, intersectionType>
        ::Quadrature(intersection, quadratureOrder, lhsBeta);
    if (type == IntegrationType::travelDistanceWeighted &&
        nOutflowIntersections > 1) {
      quadFace = SplitQuadratureRule<double>(
          quadFace,
          detail::splitPointOfInflowFaceInTriangle(intersection,
                                                   referenceBeta));
    }

    for (size_t pt=0, qsize=quadFace.size(); pt < qsize; pt++) {

      // Position of the current quadrature point in the reference element (face!)
      const FieldVector<double,dim-1>& quadFacePos = quadFace[pt].position();

      // The multiplicative factor in the integral transformation formula multiplied with outer normal
      const FieldVector<double,dim>& integrationOuterNormal =
              intersection.integrationOuterNormal(quadFacePos);

      // The multiplicative factor in the integral transformation formula -

      double integrationWeight;
      if(type == IntegrationType::normalVector ||
         type == IntegrationType::travelDistanceWeighted) {
                          // TODO: needs global geometry
        integrationWeight = detail::evaluateFactor(factor, quadFacePos)
                          * quadFace[pt].weight();
        if(type == IntegrationType::travelDistanceWeighted)
        {
          // |beta * n|*integrationweight
          integrationWeight *= fabs(lhsBeta*integrationOuterNormal);
        }
        else
          integrationWeight *= (lhsBeta*integrationOuterNormal);
      } else if(type == IntegrationType::normalSign) {
        const double integrationElement =
            intersection.geometry().integrationElement(quadFacePos);

        const FieldVector<double,dim>& centerOuterNormal =
            intersection.centerUnitOuterNormal();

        int sign = 1;
        bool signfound = false;
        for (unsigned int i=0;
           i<centerOuterNormal.size() and signfound == false;
           i++)
        {
          if (centerOuterNormal[i]<(-1e-10))
          {
            sign = -1;
            signfound = true;
          }
          else if (centerOuterNormal[i]>(1e-10))
          {
            sign = 1;
            signfound = true;
          }
        }

        // TODO: needs global geometry
        integrationWeight = sign * detail::evaluateFactor(factor, quadFacePos)
                          * quadFace[pt].weight() * integrationElement;
      }

      // position of the quadrature point within the element
      const FieldVector<double,dim> elementQuadPos =
              intersection.geometryInInside().global(quadFacePos);

      if(type == IntegrationType::travelDistanceWeighted) {
        // factor r_K(s)/|beta|
        integrationWeight *= detail::travelDistance(
            elementQuadPos,
            referenceBeta);
      }

      //////////////////////////////
      // Left Hand Side Functions //
      //////////////////////////////
      std::vector<FieldVector<double,1> > lhsValues;
      lhsLocalFiniteElement.localBasis().evaluateFunction(elementQuadPos,
                                                          lhsValues);

      ///////////////////////////////
      // Right Hand Side Functions //
      ///////////////////////////////
      std::vector<FieldVector<double,1> > rhsValues;
      rhsLocalFiniteElement.localBasis().evaluateFunction(elementQuadPos,
                                                          rhsValues);

      // Compute the actual matrix entries
      for (size_t i=0; i<nLhs; i++)
      {
        for (size_t j=0; j<nRhs; j++)
        {
          elementMatrix[i+lhsSpaceOffset][j+rhsSpaceOffset]
                  += (lhsValues[i] * rhsValues[j]) * integrationWeight;
        }
      }
    }
  }
}
};

} // end namespace detail
