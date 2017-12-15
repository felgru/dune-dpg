#include "refinedfaces.hh"

namespace Dune {
namespace detail {

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
  const int dim = Element::mydimension;
  const auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& lhsLocalFiniteElement = lhsLocalView.tree().finiteElement();
  const auto& rhsLocalFiniteElement = rhsLocalView.tree().finiteElement();

  const unsigned int nLhs(lhsLocalFiniteElement.localBasis().size());
  const unsigned int nRhs(rhsLocalFiniteElement.localBasis().size());

  typename detail::ChooseQuadrature<LhsSpace, RhsSpace, Element>::type quad
    = detail::ChooseQuadrature<LhsSpace, RhsSpace, Element>
      ::Quadrature(element, quadratureOrder);

  const auto referenceGridView =
      lhsLocalView.tree().refinedReferenceElement().leafGridView();

  const unsigned int subElementStride =
      (is_DGRefinedFiniteElement<LhsSpace>::value) ?
        lhsLocalFiniteElement.localBasis().size() : 0;

  unsigned int subElementOffset = 0;
  unsigned int subElementIndex = 0;
  for(const auto& subElement : elements(referenceGridView)) {
    const auto subGeometryInReferenceElement = subElement.geometry();
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
                  <dim, lhsType,
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
          detail::LocalFunctionEvaluation<dim, rhsType>()
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
  const int dim = Element::mydimension;
  const auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& lhsLocalFiniteElement = lhsLocalView.tree().finiteElement();
  const auto& rhsLocalFiniteElement = rhsLocalView.tree().finiteElement();

  const unsigned int nLhs(lhsLocalFiniteElement.localBasis().size());
  const unsigned int nRhs(rhsLocalFiniteElement.localBasis().size());

  const auto referenceGridView =
      lhsLocalView.tree().refinedReferenceElement().leafGridView();

  const unsigned int subElementStride =
      (is_DGRefinedFiniteElement<LhsSpace>::value) ?
        lhsLocalFiniteElement.localBasis().size() : 0;

  unsigned int subElementOffset = 0;
  unsigned int subElementIndex = 0;
  for(const auto& subElement : elements(referenceGridView))
  {
    using SubElement = std::decay_t<decltype(subElement)>;
    const auto subGeometryInReferenceElement = subElement.geometry();

    unsigned int nOutflowFaces = 0;
    for (unsigned short f = 0, fMax = subElement.subEntities(1); f < fMax; f++)
    {
      const auto face = subElement.template subEntity<1>(f);
      const FieldVector<double,dim> unitOuterNormal
          = RefinedFaceComputations<SubElement>(face, subElement, element)
              .unitOuterNormal();

      const double prod = lhsBeta * unitOuterNormal;
      if(prod > 0)
        ++nOutflowFaces;
    }

    FieldVector<double,dim> referenceBeta
        = detail::referenceBeta(geometry,
            subGeometryInReferenceElement, lhsBeta);

    for (unsigned short f = 0, fMax = subElement.subEntities(1); f < fMax; f++)
    {
      const auto face = subElement.template subEntity<1>(f);
      auto faceComputations
          = RefinedFaceComputations<SubElement>(face, subElement, element);

      using Face = std::decay_t<decltype(face)>;

      const double integrationElement = faceComputations.integrationElement();

      const FieldVector<double,dim> unitOuterNormal
          = faceComputations.unitOuterNormal();

      if(type == IntegrationType::travelDistanceWeighted &&
         lhsBeta * unitOuterNormal >= 0) {
        /* Only integrate over inflow boundaries. */
        continue;
      }

      QuadratureRule<double, 1> quadFace
        = detail::ChooseQuadrature<LhsSpace, RhsSpace, Face>
          ::Quadrature(face, quadratureOrder);
      if (type == IntegrationType::travelDistanceWeighted &&
          nOutflowFaces > 1) {
        quadFace = SplitQuadratureRule<double>(
            quadFace,
            detail::splitPointOfInflowFaceInTriangle(
                faceComputations.geometryInElement(), referenceBeta));
      }

      for (size_t pt=0, qsize=quadFace.size(); pt < qsize; pt++) {

        // Position of the current quadrature point in the reference element
        // (face!)
        const FieldVector<double,dim-1>& quadFacePos = quadFace[pt].position();

        // position of the quadrature point within the subelement
        const FieldVector<double,dim> elementQuadPosSubCell =
                faceComputations.faceToElementPosition(quadFacePos);

        // position of the quadrature point within the reference element
        const FieldVector<double,dim> elementQuadPos =
                subGeometryInReferenceElement.global(elementQuadPosSubCell);

        const FieldVector<double,dim> globalQuadPos =
                geometry.global(elementQuadPos);

        // The multiplicative factor in the integral transformation formula -

        double integrationWeight;
        if(type == IntegrationType::normalVector ||
           type == IntegrationType::travelDistanceWeighted) {
          integrationWeight = detail::evaluateFactor(factor, globalQuadPos)
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

          integrationWeight = sign
                            * detail::evaluateFactor(factor, globalQuadPos)
                            * quadFace[pt].weight() * integrationElement;
        }

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
                  <dim, EvaluationType::value,
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

}} // end namespace Dune::detail
