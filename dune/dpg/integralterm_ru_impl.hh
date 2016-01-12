template <IntegrationType type,
          class LhsSpace,
          class RhsSpace>
struct GetLocalMatrix<type, LhsSpace, RhsSpace, true, false>
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

  int level = levelOfFE<LhsSpace>::value;

  const UGGrid<dim>& referenceGrid
    = lhsLocalView.tree().refinedReferenceElement();
  auto referenceGridView = referenceGrid.leafGridView();

  assert(element.type().isTriangle() || element.type().isQuadrilateral());
  const size_t subElementStride =
    (element.type().isTriangle())
    ? lhsLocalView.globalBasis().nodeFactory().dofsPerSubTriangle
    : lhsLocalView.globalBasis().nodeFactory().dofsPerSubQuad;

  unsigned int subElementOffset = 0;
  unsigned int subElementIndex = 0;
  for(const auto& subElement : elements(referenceGridView)) {
    auto subGeometryInReferenceElement = subElement.geometry();
    for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();

      // The transposed inverse Jacobian of the map from the reference element to the element
      const auto& jacobianSub
          = subGeometryInReferenceElement.jacobianInverseTransposed(quadPos);
      const auto& jacobian = geometry.jacobianInverseTransposed
                             (subGeometryInReferenceElement.global(quadPos));

      // The multiplicative factor in the integral transformation formula
      const double integrationWeight
        = geometry.integrationElement(subGeometryInReferenceElement
                                                      .global(quadPos))
        * subGeometryInReferenceElement.integrationElement(quadPos)
        * quad[pt].weight() * detail::evaluateFactor(factor, quadPos);

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
          detail::LocalFunctionEvaluation<dim, rhsType,
                                          DomainOfIntegration::interior>()
                        (rhsLocalFiniteElement,
                         subGeometryInReferenceElement.global(quadPos),
                         geometry,
                         rhsBeta);

      // Compute the actual matrix entries
      for (unsigned int i=0; i<nLhs; i++)
      {
        for (unsigned int j=0; j<nRhs; j++)
        {
          elementMatrix[i+lhsSpaceOffset+subElementOffset]
                       [j+rhsSpaceOffset]
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

  const int nLhs(lhsLocalFiniteElement.localBasis().size());
  const int nRhs(rhsLocalFiniteElement.localBasis().size());

  const auto& gridView = lhsLocalView.globalBasis().gridView();

  int level = levelOfFE<LhsSpace>::value;

  const UGGrid<dim>& referenceGrid
    = lhsLocalView.tree().refinedReferenceElement();
  auto referenceGridView = referenceGrid.leafGridView();

  assert(element.type().isTriangle() || element.type().isQuadrilateral());
  const size_t subElementStride =
    (element.type().isTriangle())
    ? lhsLocalView.globalBasis().nodeFactory().dofsPerSubTriangle
    : lhsLocalView.globalBasis().nodeFactory().dofsPerSubQuad;

  unsigned int subElementOffset = 0;
  unsigned int subElementIndex = 0;
  for(const auto& subElement : elements(referenceGridView))
  {
    auto subGeometryInReferenceElement = subElement.geometry();

    for (auto&& intersection : intersections(gridView, subElement))
    {
      using intersectionType
        = typename std::decay<decltype(intersection)>::type;

      // compute integration element for interface
      const double integrationElement
        = (   geometry.global(intersection.geometry().global({1}))
            - geometry.global(intersection.geometry().global({0})) )
          .two_norm();

      // TODO: Do we really want to have a transport quadrature rule
      //       on the faces, if one of the FE spaces is a transport space?
      typename detail::ChooseQuadrature<LhsSpace,
                                        RhsSpace,
                                        intersectionType>::type quadFace
        = detail::ChooseQuadrature<LhsSpace, RhsSpace, intersectionType>
          ::Quadrature(intersection, quadratureOrder, lhsBeta);

      for (size_t pt=0, qsize=quadFace.size(); pt < qsize; pt++) {

        // Position of the current quadrature point in the reference element (face!)
        const FieldVector<double,dim-1>& quadFacePos = quadFace[pt].position();
        const FieldVector<double,dim>& quadFacePosInReferenceElement
          = intersection.geometry().global(quadFacePos);

        // The multiplicative factor in the integral transformation formula multiplied with outer normal
        // TODO: needs global geometry
        const FieldVector<double,dim>& integrationOuterNormal =
              // TODO: This is missing the reference → physical domain transformation (is missing factor 1/n; do we apply this later?):
              intersection.integrationOuterNormal(quadFacePos);

        // The multiplicative factor in the integral transformation formula -

        double integrationWeight;
        if(type == IntegrationType::normalVector) {
          integrationWeight = (lhsBeta*integrationOuterNormal)
                            // TODO: needs global geometry
                            * detail::evaluateFactor(factor, quadFacePos)
                            * quadFace[pt].weight()
                            * integrationElement;
        } else if(type == IntegrationType::normalSign) {
          const FieldVector<double,dim>& centerOuterNormal =
              // TODO: needs global geometry
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

        // position of the quadrature point within the subelement
        const FieldVector<double,dim> elementQuadPosSubCell =
                intersection.geometryInInside().global(quadFacePos);

        // position of the quadrature point within the reference element
        const FieldVector<double,dim> elementQuadPos =
                subGeometryInReferenceElement.global(elementQuadPosSubCell);


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
        std::vector<FieldVector<double,1> > rhsValues;
        rhsLocalFiniteElement.localBasis()
            .evaluateFunction(elementQuadPos, rhsValues);

        // Compute the actual matrix entries
        for (size_t i=0; i<nLhs; i++)
        {
          for (size_t j=0; j<nRhs; j++)
          {
            elementMatrix[i+lhsSpaceOffset+subElementOffset]
                         [j+rhsSpaceOffset]
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
