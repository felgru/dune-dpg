namespace detail {

template <IntegrationType type,
          class LhsSpace,
          class RhsSpace>
struct GetLocalMatrix<type, LhsSpace, RhsSpace, true, true>
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

  const auto& referenceGrid
    = lhsLocalView.tree().refinedReferenceElement();
  auto referenceGridView = referenceGrid.leafGridView();

  const size_t subElementStride =
      (is_DGRefinedFiniteElement<LhsSpace>::value) ?
        lhsLocalFiniteElement.localBasis().size() : 0;

  unsigned int subElementOffset = 0;
  unsigned int subElementIndex = 0;
  for(const auto& subElement : elements(referenceGridView)) {
    auto subGeometryInReferenceElement = subElement.geometry();
    for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();
     // Global position of the current quadrature point
     const FieldVector<double,dim>& globalQuadPos
          = geometry.global(subGeometryInReferenceElement.global(quadPos));

      // The multiplicative factor in the integral transformation formula
      const double integrationWeight
        = geometry.integrationElement(subGeometryInReferenceElement
                                                      .global(quadPos))
        * subGeometryInReferenceElement.integrationElement(quadPos)
        * quad[pt].weight() * detail::evaluateFactor(factor, globalQuadPos);

      //////////////////////////////
      // Left hand side Functions //
      //////////////////////////////
      constexpr auto lhsType = (type == IntegrationType::valueValue ||
                                type == IntegrationType::valueGrad)
                               ? EvaluationType::value : EvaluationType::grad;

      std::vector<FieldVector<double,1> > lhsValues =
          detail::LocalRefinedFunctionEvaluation
                  <dim, lhsType, DomainOfIntegration::interior,
                   is_ContinuouslyRefinedFiniteElement<LhsSpace>::value>()
                        (lhsLocalFiniteElement,
                         subElementIndex,
                         quadPos,
                         geometry,
                         subGeometryInReferenceElement,
                         lhsBeta);

      ///////////////////////////////
      // Right hand side Functions //
      ///////////////////////////////
      constexpr auto rhsType = (type == IntegrationType::valueValue ||
                                type == IntegrationType::gradValue)
                               ? EvaluationType::value : EvaluationType::grad;

      std::vector<FieldVector<double,1> > rhsValues =
          detail::LocalRefinedFunctionEvaluation
                  <dim, rhsType, DomainOfIntegration::interior,
                   is_ContinuouslyRefinedFiniteElement<RhsSpace>::value>()
                        (rhsLocalFiniteElement,
                         subElementIndex,
                         quadPos,
                         geometry,
                         subGeometryInReferenceElement,
                         rhsBeta);

      // Compute the actual matrix entries
      for (unsigned int i=0; i<nLhs; i++)
      {
        for (unsigned int j=0; j<nRhs; j++)
        {
          elementMatrix[i+lhsSpaceOffset+subElementOffset]
                       [j+rhsSpaceOffset+subElementOffset]
                  += (lhsValues[i] * rhsValues[j]) * integrationWeight;
        }
      }
    }
    if(is_DGRefinedFiniteElement<LhsSpace>::value)
      subElementOffset += subElementStride;
    subElementIndex++;
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
  auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& lhsLocalFiniteElement = lhsLocalView.tree().finiteElement();
  const auto& rhsLocalFiniteElement = rhsLocalView.tree().finiteElement();

  const unsigned int nLhs(lhsLocalFiniteElement.localBasis().size());
  const unsigned int nRhs(rhsLocalFiniteElement.localBasis().size());

  const auto& gridView = lhsLocalView.globalBasis().gridView();

  const auto& referenceGrid
    = lhsLocalView.tree().refinedReferenceElement();
  auto referenceGridView = referenceGrid.leafGridView();

  const size_t subElementStride =
      (is_DGRefinedFiniteElement<LhsSpace>::value) ?
        lhsLocalFiniteElement.localBasis().size() : 0;

  unsigned int subElementOffset = 0;
  unsigned int subElementIndex = 0;
  for(const auto& subElement : elements(referenceGridView))
  {
    auto subGeometryInReferenceElement = subElement.geometry();

    unsigned int nInflowIntersections = 0;
    unsigned int nOutflowIntersections = 0;
    for (auto&& intersection : intersections(referenceGridView, subElement))
    {
      using intersectionType
        = typename std::decay<decltype(intersection)>::type;

      const FieldVector<double,dim> globalCorner0
        = geometry.global(intersection.geometry().global({0}));
      const FieldVector<double,dim> globalCorner1
        = geometry.global(intersection.geometry().global({1}));
      // compute integration element for interface
      const double integrationElement
        = (globalCorner1 - globalCorner0).two_norm();

      static_assert(dim==2, "Computation of unit outer normal for subcell"
                            " only implemented in 2d!");
      /* This won't work for curvilinear elements, but they don't seem
       * to be supported by UG anyway. */
      FieldVector<double,dim> unitOuterNormal
        = { (globalCorner1[1] - globalCorner0[1])
          , (globalCorner0[0] - globalCorner1[0]) };
      unitOuterNormal /= unitOuterNormal.two_norm();

      double prod = lhsBeta * unitOuterNormal;
      if(prod > 0)
        ++nOutflowIntersections;
      else if (prod < 0)
        ++nInflowIntersections;
    }

    FieldVector<double,dim> referenceBeta
        = detail::referenceBeta(geometry,
            subGeometryInReferenceElement, lhsBeta);

    for (auto&& intersection : intersections(referenceGridView, subElement))
    {
      using intersectionType
        = typename std::decay<decltype(intersection)>::type;

      const FieldVector<double,dim> globalCorner0
        = geometry.global(intersection.geometry().global({0}));
      const FieldVector<double,dim> globalCorner1
        = geometry.global(intersection.geometry().global({1}));
      // compute integration element for interface
      const double integrationElement
        = (globalCorner1 - globalCorner0).two_norm();

      static_assert(dim==2, "Computation of unit outer normal for subcell"
                            " only implemented in 2d!");
      /* This won't work for curvilinear elements, but they don't seem
       * to be supported by UG anyway. */
      FieldVector<double,dim> unitOuterNormal
        = { (globalCorner1[1] - globalCorner0[1])
          , (globalCorner0[0] - globalCorner1[0]) };
      unitOuterNormal /= unitOuterNormal.two_norm();

      if(type == IntegrationType::travelDistanceWeighted &&
         lhsBeta * unitOuterNormal >= 0) {
        /* Only integrate over inflow boundaries. */
        continue;
      }

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
        // const FieldVector<double,dim>& quadFacePosInReferenceElement
        //   = intersection.geometry().global(quadFacePos);

        // The multiplicative factor in the integral transformation formula -

        double integrationWeight;
        if(type == IntegrationType::normalVector ||
           type == IntegrationType::travelDistanceWeighted) {
                            // TODO: needs global geometry
          integrationWeight = detail::evaluateFactor(factor, quadFacePos)
                            * quadFace[pt].weight()
                            * integrationElement;
          // TODO: scale lhsBeta to length 1
          if(type == IntegrationType::travelDistanceWeighted)
            integrationWeight *= fabs(lhsBeta*unitOuterNormal);
          else
            integrationWeight *= (lhsBeta*unitOuterNormal);
        } else if(type == IntegrationType::normalSign) {
          int sign = 1;
          bool signfound = false;
          for (unsigned int i=0;
             i < unitOuterNormal.size() and signfound == false;
             i++)
          {
            if (unitOuterNormal[i]<(-1e-10))
            {
              sign = -1;
              signfound = true;
            }
            else if (unitOuterNormal[i]>(1e-10))
            {
              sign = 1;
              signfound = true;
            }
          }

          // TODO: needs global geometry
          integrationWeight = sign * detail::evaluateFactor(factor, quadFacePos)
                            * quadFace[pt].weight() * integrationElement;
        }

        // position of the quadrature point within the subelement
        const FieldVector<double,dim> elementQuadPosSubCell =
                intersection.geometryInInside().global(quadFacePos);

        if(type == IntegrationType::travelDistanceWeighted) {
          integrationWeight *= detail::travelDistance(
              elementQuadPosSubCell,
              referenceBeta);
        }

        //////////////////////////////
        // Left Hand Side Functions //
        //////////////////////////////
        std::vector<FieldVector<double,1> > lhsValues =
          detail::LocalRefinedFunctionEvaluation
                  <dim, EvaluationType::value, DomainOfIntegration::interior,
                   is_ContinuouslyRefinedFiniteElement<LhsSpace>::value>()
                        (lhsLocalFiniteElement,
                         subElementIndex,
                         elementQuadPosSubCell,
                         geometry,
                         subGeometryInReferenceElement,
                         lhsBeta);

        ///////////////////////////////
        // Right Hand Side Functions //
        ///////////////////////////////
        std::vector<FieldVector<double,1> > rhsValues =
          detail::LocalRefinedFunctionEvaluation
                  <dim, EvaluationType::value, DomainOfIntegration::interior,
                   is_ContinuouslyRefinedFiniteElement<RhsSpace>::value>()
                        (rhsLocalFiniteElement,
                         subElementIndex,
                         elementQuadPosSubCell,
                         geometry,
                         subGeometryInReferenceElement,
                         rhsBeta);

        // Compute the actual matrix entries
        for (size_t i=0; i<nLhs; i++)
        {
          for (size_t j=0; j<nRhs; j++)
          {
            elementMatrix[i+lhsSpaceOffset+subElementOffset]
                         [j+rhsSpaceOffset+subElementOffset]
                    += (lhsValues[i] * rhsValues[j]) * integrationWeight;
          }
        }
      }
    }
    if(is_DGRefinedFiniteElement<LhsSpace>::value)
      subElementOffset += subElementStride;
    subElementIndex++;
  }
}
};

} // end namespace detail
