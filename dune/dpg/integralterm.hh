// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_INTEGRALTERM_HH
#define DUNE_DPG_INTEGRALTERM_HH

#include <tuple>
#include <vector>
#include <type_traits>

#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/istl/matrix.hh>

#include <dune/common/std/final.hh>

#include <dune/functions/functionspacebases/interpolate.hh>

#include "assemble_types.hh"
#include "traits.hh"

namespace Dune {

  /**
   * class IntegralTerm
   *
   */
  template <IntegrationType type,
            DomainOfIntegration domain_of_integration =
                DomainOfIntegration::interior>
  class IntegralTerm
  {
  public:

    /* TODO: make this a template argument. */
    static const int dim = 2;

    IntegralTerm () = delete;

    IntegralTerm (double constant_factor = 1,
                  FieldVector<double, dim> lhsBeta = {1,1},
                  FieldVector<double, dim> rhsBeta = {1,1})
        : constant_factor(constant_factor),
          lhsBeta(lhsBeta),
          rhsBeta(rhsBeta)
    {};

    /** Compute the stiffness matrix for a single element */
    template <class LhsLocalView,
              class RhsLocalView,
              class MatrixType>
    void getLocalMatrix(const LhsLocalView& lhsLocalView,
                        const RhsLocalView& rhsLocalView,
                        MatrixType& elementMatrix,
                        size_t lhsSpaceOffset,
                        size_t rhsSpaceOffset) const;

  private:

    double constant_factor;
    FieldVector<double, dim> lhsBeta;
    FieldVector<double, dim> rhsBeta;

  };

template<size_t lhsSpaceIndex,
         size_t rhsSpaceIndex,
         IntegrationType integrationType,
         DomainOfIntegration domainOfIntegration,
         EnableIf<std::integral_constant<bool,
                     integrationType == IntegrationType::valueValue
                  || integrationType == IntegrationType::normalSign> >...
        >
auto make_IntegralTerm(double c)
    -> std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                  std::integral_constant<size_t, rhsSpaceIndex>,
                  IntegralTerm<integrationType, domainOfIntegration> >
{
  return std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                std::integral_constant<size_t, rhsSpaceIndex>,
                IntegralTerm<integrationType, domainOfIntegration> >
         ({},{},
          IntegralTerm<integrationType, domainOfIntegration>(c));
}

template<size_t lhsSpaceIndex,
         size_t rhsSpaceIndex,
         IntegrationType integrationType,
         DomainOfIntegration domainOfIntegration,
         EnableIf<std::integral_constant<bool,
                     integrationType == IntegrationType::gradValue
                  || integrationType == IntegrationType::valueGrad
                  || integrationType == IntegrationType::gradGrad
                  || integrationType == IntegrationType::normalVector> >...
        >
auto make_IntegralTerm(double c, FieldVector<double, 2> beta)
    -> std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                  std::integral_constant<size_t, rhsSpaceIndex>,
                  IntegralTerm<integrationType, domainOfIntegration> >
{
  return std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                std::integral_constant<size_t, rhsSpaceIndex>,
                IntegralTerm<integrationType, domainOfIntegration> >
         ({},{},
          IntegralTerm<integrationType, domainOfIntegration>(c, beta, beta));
}

template<size_t lhsSpaceIndex,
         size_t rhsSpaceIndex,
         IntegrationType integrationType,
         DomainOfIntegration domainOfIntegration,
         EnableIf<std::integral_constant<bool,
                     integrationType == IntegrationType::gradGrad> >...
        >
auto make_IntegralTerm(double c,
                       FieldVector<double, 2> lhsBeta,
                       FieldVector<double, 2> rhsBeta)
    -> std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                  std::integral_constant<size_t, rhsSpaceIndex>,
                  IntegralTerm<integrationType, domainOfIntegration> >
{
  return std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                std::integral_constant<size_t, rhsSpaceIndex>,
                IntegralTerm<integrationType, domainOfIntegration> >
         ({},{},
          IntegralTerm<integrationType, domainOfIntegration>(c, lhsBeta,
                                                                rhsBeta));
}

template<IntegrationType type, DomainOfIntegration domain_of_integration>
template <class LhsLocalView,
          class RhsLocalView,
          class MatrixType>
void IntegralTerm<type, domain_of_integration>
     ::getLocalMatrix(
        const LhsLocalView& lhsLocalView,
        const RhsLocalView& rhsLocalView,
        MatrixType& elementMatrix,
        size_t lhsSpaceOffset,
        size_t rhsSpaceOffset) const
{
  // Get the grid element from the local FE basis view
  using Element = typename std::remove_pointer<LhsLocalView>::type::Element;
  const Element& element = lhsLocalView->element();

  const int dim = Element::dimension;
  auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& lhsLocalFiniteElement = lhsLocalView->tree().finiteElement();
  const auto& rhsLocalFiniteElement = rhsLocalView->tree().finiteElement();

  const int nLhs(lhsLocalFiniteElement.localBasis().size());
  const int nRhs(rhsLocalFiniteElement.localBasis().size());

  // Order for the quadrature rule
  /* TODO: can probably be one less for gradients. */
  int order = 2*(dim*lhsLocalFiniteElement.localBasis().order()-1);

  ////////////////////////////
  // Assemble interior terms
  ////////////////////////////
  if(domain_of_integration == DomainOfIntegration::interior) {

  // Get a quadrature rule
  const QuadratureRule<double, dim>& quad =
          QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (size_t pt=0; pt < quad.size(); pt++) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The transposed inverse Jacobian of the map from the reference element to the element
    const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPos);

    //////////////////////////////
    // Left hand side Functions //
    //////////////////////////////
    // values of the shape functions
    std::vector<FieldVector<double,1> > lhsValues;
    lhsLocalFiniteElement.localBasis().evaluateFunction(quadPos, lhsValues);

    /* TODO: only compute what is necessary */
    // The gradients of the shape functions on the reference element
    std::vector<FieldMatrix<double,1,dim> > lhsReferenceGradients;
    lhsLocalFiniteElement.localBasis()
            .evaluateJacobian(quadPos, lhsReferenceGradients);

    // Compute the shape function gradients on the real element
    std::vector<FieldVector<double,dim> >
            lhsGradients(lhsReferenceGradients.size());
        for (size_t i=0; i<lhsGradients.size(); i++)
      jacobian.mv(lhsReferenceGradients[i][0], lhsGradients[i]);

    ///////////////////////////////
    // Right hand side Functions //
    ///////////////////////////////
    // values of the shape functions
    std::vector<FieldVector<double,1> > rhsValues;
    rhsLocalFiniteElement.localBasis()
            .evaluateFunction(quadPos, rhsValues);

    // The gradients of the shape functions on the reference element
    std::vector<FieldMatrix<double,1,dim> > rhsReferenceGradients;
    rhsLocalFiniteElement.localBasis()
            .evaluateJacobian(quadPos, rhsReferenceGradients);

    // Compute the shape function gradients on the real element
    std::vector<FieldVector<double,dim> >
            rhsGradients(rhsReferenceGradients.size());
    for (size_t i=0; i<rhsGradients.size(); i++)
      jacobian.mv(rhsReferenceGradients[i][0], rhsGradients[i]);

    // Compute the actual matrix entries
    for (size_t i=0; i<nLhs; i++)
    {
      for (size_t j=0; j<nRhs; j++)
      {
        static_assert(type == IntegrationType::valueValue
                   || type == IntegrationType::gradValue
                   || type == IntegrationType::valueGrad
                   || type == IntegrationType::gradGrad
                   || type == IntegrationType::normalVector
                   || type == IntegrationType::normalSign,
                   "Use of unknown IntegrationType.");
        static_assert(domain_of_integration != DomainOfIntegration::interior
                      || type == IntegrationType::valueValue
                      || type == IntegrationType::gradValue
                      || type == IntegrationType::valueGrad
                      || type == IntegrationType::gradGrad,
                      "IntegrationType not implemented on interior.");
        if(type == IntegrationType::valueValue) {
        elementMatrix[i+lhsSpaceOffset][j+rhsSpaceOffset]
                += (lhsValues[i] * rhsValues[j]) * constant_factor
                   * quad[pt].weight() * integrationElement;
        } else if(type == IntegrationType::valueGrad) {
        elementMatrix[i+lhsSpaceOffset][j+rhsSpaceOffset]
                += lhsValues[i] * (rhsBeta*rhsGradients[j]) * constant_factor
                   * quad[pt].weight() * integrationElement;
        } else if(type == IntegrationType::gradValue) {
        elementMatrix[i+lhsSpaceOffset][j+rhsSpaceOffset]
                += (lhsBeta*lhsGradients[i]) * rhsValues[j] * constant_factor
                   * quad[pt].weight() * integrationElement;
        } else if(type == IntegrationType::gradGrad) {
        elementMatrix[i+lhsSpaceOffset][j+rhsSpaceOffset]
                += (lhsBeta*lhsGradients[i]) * (rhsBeta*rhsGradients[j])
                   * constant_factor * quad[pt].weight() * integrationElement;
        }
      }
    }
  }

  } else {
  ////////////////////////////
  // Assemble boundary terms
  ////////////////////////////

  const auto& gridView = lhsLocalView->globalBasis().gridView();

  for (auto&& intersection : intersections(gridView, element))
  {
    const QuadratureRule<double, dim-1>& quadFace =
            QuadratureRules<double, dim-1>::rule(intersection.type(), order);
    // Loop over all quadrature points
    for (size_t pt=0; pt < quadFace.size(); pt++) {

    // Position of the current quadrature point in the reference element (face!)
    const FieldVector<double,dim-1>& quadFacePos = quadFace[pt].position();

    // The multiplicative factor in the integral transformation formula multiplied with outer normal
    const FieldVector<double,dim>& integrationOuterNormal =
            intersection.integrationOuterNormal(quadFacePos);

    // The multiplicative factor in the integral transformation formula -
    // TODO this way, the computation is absolutely inefficient!
    const double integrationElement = integrationOuterNormal.two_norm();
   //const double integrationElement = intersection.geometry().integrationElement(quadFacePos); TODO ???

   const FieldVector<double,dim>& centerOuterNormal =
            intersection.centerUnitOuterNormal();

    int sign = 1;
    bool signfound = false;
    for (unsigned int i=0;
         i<centerOuterNormal.size() and signfound == false;
         i++)
    {
      if (centerOuterNormal[i]<0)
      {
        sign = -1;
        signfound = true;
      }
      else if (centerOuterNormal[i]>0)
      {
        sign = 1;
        signfound = true;
      }
    }


                // position of the quadrature point within the element
    const FieldVector<double,dim> elementQuadPos =
            intersection.geometryInInside().global(quadFacePos);


    //////////////////////////////
    // Left Hand Side Functions //
    //////////////////////////////
    // values of the shape functions
    std::vector<FieldVector<double,1> > lhsValues;
    lhsLocalFiniteElement.localBasis().evaluateFunction(elementQuadPos,
                                                        lhsValues);

    ///////////////////////////////
    // Right Hand Side Functions //
    ///////////////////////////////
    // values of the shape functions
    std::vector<FieldVector<double,1> > rhsValues;
    rhsLocalFiniteElement.localBasis().evaluateFunction(elementQuadPos,
                                                        rhsValues);

    // Compute the actual matrix entries
    for (size_t i=0; i<nLhs; i++)
    {
      for (size_t j=0; j<nRhs; j++)
      {
        static_assert(type == IntegrationType::valueValue
                   || type == IntegrationType::gradValue
                   || type == IntegrationType::valueGrad
                   || type == IntegrationType::gradGrad
                   || type == IntegrationType::normalVector
                   || type == IntegrationType::normalSign,
                   "Use of unknown IntegrationType.");
        static_assert(domain_of_integration != DomainOfIntegration::face
                      || type == IntegrationType::normalVector
                      || type == IntegrationType::normalSign,
                      "IntegrationType not implemented on boundary.");
        if(type == IntegrationType::normalVector) {
        elementMatrix[i+lhsSpaceOffset][j+rhsSpaceOffset]
                += ((lhsBeta*integrationOuterNormal) * constant_factor
                    * lhsValues[i] * rhsValues[j]) * quadFace[pt].weight();
        } else if(type == IntegrationType::normalSign) {
        elementMatrix[i+lhsSpaceOffset][j+rhsSpaceOffset]
                += (sign * constant_factor * lhsValues[i] * rhsValues[j])
                    * quadFace[pt].weight() * integrationElement;
        }
      }
    }
    }
  }
  }
}

} // end namespace Dune

#endif // DUNE_DPG_INTEGRALTERM_HH
