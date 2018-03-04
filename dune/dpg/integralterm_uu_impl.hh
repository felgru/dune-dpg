#include "faces.hh"

namespace Dune {
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
          class LocalFactor,
          class DirectionType>
inline static void interiorImpl(const LhsLocalView& lhsLocalView,
                                const RhsLocalView& rhsLocalView,
                                MatrixType& elementMatrix,
                                size_t lhsSpaceOffset,
                                size_t rhsSpaceOffset,
                                unsigned int quadratureOrder,
                                const Element& element,
                                const LocalFactor& localFactor,
                                const DirectionType& lhsBeta,
                                const DirectionType& rhsBeta)
{
  constexpr int dim = Element::mydimension;
  const auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& lhsLocalFiniteElement = lhsLocalView.tree().finiteElement();
  const auto& rhsLocalFiniteElement = rhsLocalView.tree().finiteElement();

  const unsigned int nLhs(lhsLocalFiniteElement.size());
  const unsigned int nRhs(rhsLocalFiniteElement.size());

  typename detail::ChooseQuadrature<LhsSpace, RhsSpace, Element>::type quad
    = detail::ChooseQuadrature<LhsSpace, RhsSpace, Element>
      ::Quadrature(element, quadratureOrder);

  for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();
    // Global position of the current quadrature point

    // The multiplicative factor in the integral transformation formula
    const double integrationWeight = geometry.integrationElement(quadPos)
                                   * quad[pt].weight()
                                   * localFactor(quadPos);


    //////////////////////////////
    // Left hand side Functions //
    //////////////////////////////
    constexpr auto lhsType = (type == IntegrationType::valueValue ||
                              type == IntegrationType::valueGrad)
                             ? EvaluationType::value : EvaluationType::grad;

    const std::vector<FieldVector<double,1> > lhsValues =
        detail::LocalFunctionEvaluation<dim, lhsType>()
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

    const std::vector<FieldVector<double,1> > rhsValues =
        detail::LocalFunctionEvaluation<dim, rhsType>()
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
          class LocalFactor,
          class DirectionType>
inline static void
faceImpl(const LhsLocalView& lhsLocalView,
         const RhsLocalView& rhsLocalView,
         MatrixType& elementMatrix,
         size_t lhsSpaceOffset,
         size_t rhsSpaceOffset,
         unsigned int quadratureOrder,
         const Element& element,
         const LocalFactor& localFactor,
         const DirectionType& lhsBeta,
         const DirectionType& /* rhsBeta */)
{
  constexpr int dim = Element::mydimension;

  // Get set of shape functions for this element
  const auto& lhsLocalFiniteElement = lhsLocalView.tree().finiteElement();
  const auto& rhsLocalFiniteElement = rhsLocalView.tree().finiteElement();

  const unsigned int nLhs(lhsLocalFiniteElement.size());
  const unsigned int nRhs(rhsLocalFiniteElement.size());

  unsigned int nOutflowFaces = 0;
  for (unsigned short f = 0, fMax = element.subEntities(1); f < fMax; f++)
  {
    auto face = element.template subEntity<1>(f);
    const double prod = lhsBeta
      * FaceComputations<Element>(face, element).unitOuterNormal();
    if(prod > 0)
      ++nOutflowFaces;
  }

  const auto geometry = element.geometry();

  FieldVector<double,dim> referenceBeta;
  {
    const auto& jacobianInverse = geometry.jacobianInverseTransposed({0., 0.});
    jacobianInverse.mtv(lhsBeta, referenceBeta);
  }

  for (unsigned short f = 0, fMax = element.subEntities(1); f < fMax; f++)
  {
    auto face = element.template subEntity<1>(f);
    auto faceComputations = FaceComputations<Element>(face, element);
    if(type == IntegrationType::travelDistanceWeighted &&
       lhsBeta * faceComputations.unitOuterNormal() >= 0) {
      /* Only integrate over inflow boundaries. */
      continue;
    }

    using Face = std::decay_t<decltype(face)>;
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

      // position of the quadrature point within the element
      const FieldVector<double,dim> elementQuadPos =
              faceComputations.faceToElementPosition(quadFacePos);

      // The multiplicative factor in the integral transformation
      // formula multiplied with outer normal
      const FieldVector<double,dim> integrationOuterNormal =
              faceComputations.integrationOuterNormal();

      // The multiplicative factor in the integral transformation formula
      double integrationWeight;
      if(type == IntegrationType::normalVector ||
         type == IntegrationType::travelDistanceWeighted) {
        integrationWeight = localFactor(elementQuadPos)
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
            face.geometry().integrationElement(quadFacePos);

        const FieldVector<double,dim>& centerOuterNormal =
               faceComputations.unitOuterNormal();

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

        integrationWeight = sign
                          * localFactor(elementQuadPos)
                          * quadFace[pt].weight() * integrationElement;
      }

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

}} // end namespace Dune::detail
