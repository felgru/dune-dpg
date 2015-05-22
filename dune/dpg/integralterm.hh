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
   * \brief This class describes an integral term.
   *
   * This is the essential building block from which BilinearForm and
   * InnerProduct are built.
   *
   * \tparam integrationType  the form of the integrand
   * \tparam domainOfIntegration
   * \tparam FactorType     the type of the factor with which
   *                        we multiply the integrand
   * \tparam DirectionType  the type of the transport directions
   */
  template <IntegrationType type,
            DomainOfIntegration domain_of_integration,
            class FactorType,
            class DirectionType = FieldVector<double, 2> >
  class IntegralTerm
  {
  public:

    IntegralTerm () = delete;

    /**
     * \brief constructor for IntegralTerm
     *
     * \note For your convenience, use make_IntegralTerm() instead.
     */
    IntegralTerm (FactorType factor = 1,
                  DirectionType lhsBeta = {1,1},
                  DirectionType rhsBeta = {1,1})
        : factor(factor),
          lhsBeta(lhsBeta),
          rhsBeta(rhsBeta)
    {};

    /**
     * \brief Compute the stiffness matrix for a single element.
     *
     * The local integrals will be added with the given offsets
     * to \p elementMatrix.
     *
     * \pre The localViews have to be bound to the same element.
     *
     * \param[in]     lhsLocalView    local view of the left space
     * \param[in]     lhsLocalView    local view of the right space
     * \param[in,out] elementMatrix   the local system matrix
     * \param         lhsSpaceOffset  row offset for the left space
     * \param         rhsSpaceOffset  column offset for the right space
     */
    template <class LhsLocalView,
              class RhsLocalView,
              class MatrixType>
    void getLocalMatrix(const LhsLocalView& lhsLocalView,
                        const RhsLocalView& rhsLocalView,
                        MatrixType& elementMatrix,
                        size_t lhsSpaceOffset,
                        size_t rhsSpaceOffset) const;

  private:

    FactorType factor;
    DirectionType lhsBeta;
    DirectionType rhsBeta;

  };

/**
 * \brief Creates a Tuple of an IntegralTerm and the indices
 *        of both spaces involved.
 *
 * \param c  the factor with which we multiply the integrand
 * \tparam lhsSpaceIndex the index of the left space
 * \tparam rhsSpaceIndex the index of the right space
 * \tparam integrationType  the form of the integrand
 * \tparam domainOfIntegration
 * \tparam FactorType  the type of the factor \p c
 */
template<size_t lhsSpaceIndex,
         size_t rhsSpaceIndex,
         IntegrationType integrationType,
         DomainOfIntegration domainOfIntegration,
         class FactorType,
         EnableIf<std::integral_constant<bool,
                     integrationType == IntegrationType::valueValue
                  || integrationType == IntegrationType::normalSign> >...
        >
auto make_IntegralTerm(FactorType c)
    -> std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                  std::integral_constant<size_t, rhsSpaceIndex>,
                  IntegralTerm<integrationType, domainOfIntegration,
                               FactorType> >
{
  return std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                std::integral_constant<size_t, rhsSpaceIndex>,
                IntegralTerm<integrationType, domainOfIntegration,
                             FactorType> >
         ({},{},
          IntegralTerm<integrationType, domainOfIntegration, FactorType>(c));
}

/**
 * \brief Creates a Tuple of an IntegralTerm and the indices
 *        of both spaces involved.
 *
 * \param c     the factor with which we multiply the integrand
 * \param beta  the transport direction
 * \tparam lhsSpaceIndex the index of the left space
 * \tparam rhsSpaceIndex the index of the right space
 * \tparam integrationType  the form of the integrand
 * \tparam domainOfIntegration
 * \tparam FactorType     the type of the factor \p c
 * \tparam DirectionType  the type of the transport direction \p beta
 */
template<size_t lhsSpaceIndex,
         size_t rhsSpaceIndex,
         IntegrationType integrationType,
         DomainOfIntegration domainOfIntegration,
         class FactorType, class DirectionType,
         EnableIf<std::integral_constant<bool,
                     integrationType == IntegrationType::gradValue
                  || integrationType == IntegrationType::valueGrad
                  || integrationType == IntegrationType::gradGrad
                  || integrationType == IntegrationType::normalVector> >...
        >
auto make_IntegralTerm(FactorType c, DirectionType beta)
    -> std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                  std::integral_constant<size_t, rhsSpaceIndex>,
                  IntegralTerm<integrationType, domainOfIntegration,
                               FactorType, DirectionType> >
{
  return std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                std::integral_constant<size_t, rhsSpaceIndex>,
                IntegralTerm<integrationType, domainOfIntegration,
                             FactorType, DirectionType> >
         ({},{},
          IntegralTerm<integrationType, domainOfIntegration,
                       FactorType, DirectionType>(c, beta, beta));
}

/**
 * \brief Creates a Tuple of an IntegralTerm and the indices
 *        of both spaces involved.
 *
 * \param c        the factor with which we multiply the integrand
 * \param lhsbeta  the transport direction for the left space
 * \param rhsbeta  the transport direction for the right space
 * \tparam lhsSpaceIndex the index of the left space
 * \tparam rhsSpaceIndex the index of the right space
 * \tparam integrationType  the form of the integrand
 * \tparam domainOfIntegration
 * \tparam FactorType     the type of the factor \p c
 * \tparam DirectionType  the type of the transport directions
 */
template<size_t lhsSpaceIndex,
         size_t rhsSpaceIndex,
         IntegrationType integrationType,
         DomainOfIntegration domainOfIntegration,
         class FactorType, class DirectionType,
         EnableIf<std::integral_constant<bool,
                     integrationType == IntegrationType::gradGrad> >...
        >
auto make_IntegralTerm(FactorType c,
                       DirectionType lhsBeta,
                       DirectionType rhsBeta)
    -> std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                  std::integral_constant<size_t, rhsSpaceIndex>,
                  IntegralTerm<integrationType, domainOfIntegration,
                               FactorType, DirectionType> >
{
  return std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                std::integral_constant<size_t, rhsSpaceIndex>,
                IntegralTerm<integrationType, domainOfIntegration,
                             FactorType, DirectionType> >
         ({},{},
          IntegralTerm<integrationType, domainOfIntegration,
                       FactorType, DirectionType>(c, lhsBeta, rhsBeta));
}

namespace detail {
    /* We need to make this a class, as partial specializations of
     * function templates are not allowed. */
template<int dim, EvaluationType type,
         DomainOfIntegration domain_of_integration>
struct LocalFunctionEvaluation {

  template <class LocalFiniteElement, class Geometry>
  std::vector<FieldVector<double,1> >
  operator() (const LocalFiniteElement& localFiniteElement,
              const FieldVector<double, dim>& quadPos,
              const Geometry& geometry,
              const FieldVector<double, dim>& beta) const;
};

template<int dim, DomainOfIntegration domain_of_integration>
struct LocalFunctionEvaluation<dim, EvaluationType::value,
                               domain_of_integration> {

  template <class LocalFiniteElement, class Geometry>
  std::vector<FieldVector<double,1> > operator()
                      (const LocalFiniteElement& localFiniteElement,
                       const FieldVector<double, dim>& quadPos,
                       const Geometry& geometry,
                       const FieldVector<double, dim>&) const
  {
    // values of the shape functions
    std::vector<FieldVector<double,1> > values;
    localFiniteElement.localBasis().evaluateFunction(quadPos, values);
    return values;
  }
};

template<int dim, DomainOfIntegration domain_of_integration>
struct LocalFunctionEvaluation<dim, EvaluationType::grad,
                               domain_of_integration> {

  template <class LocalFiniteElement, class Geometry>
  std::vector<FieldVector<double,1> > operator()
                      (const LocalFiniteElement& localFiniteElement,
                       const FieldVector<double, dim> & quadPos,
                       const Geometry& geometry,
                       const FieldVector<double, dim>& beta) const
  {
    const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);
    // The gradients of the shape functions on the reference element
    std::vector<FieldMatrix<double,1,dim> > referenceGradients;
    localFiniteElement.localBasis()
            .evaluateJacobian(quadPos, referenceGradients);

    // Compute the shape function gradients on the real element
    std::vector<FieldVector<double, 1> >
            derivatives(referenceGradients.size());
    for (size_t i=0, i_max=referenceGradients.size(); i<i_max; i++)
    {
      FieldVector<double,dim> gradient;
      jacobian.mv(referenceGradients[i][0], gradient);
      derivatives[i] = beta * gradient;
    }

    return derivatives;
  }
};

}

template<IntegrationType type, DomainOfIntegration domain_of_integration,
         class FactorType, class DirectionType>
template <class LhsLocalView,
          class RhsLocalView,
          class MatrixType>
void IntegralTerm<type, domain_of_integration, FactorType, DirectionType>
     ::getLocalMatrix(
        const LhsLocalView& lhsLocalView,
        const RhsLocalView& rhsLocalView,
        MatrixType& elementMatrix,
        size_t lhsSpaceOffset,
        size_t rhsSpaceOffset) const
{
  static_assert(std::is_arithmetic<typename std::decay<FactorType>::type
                                  >::value
             && std::is_same<typename std::decay<DirectionType>::type,
                             FieldVector<double, 2>
                            >::value,
             "getLocalMatrix only implemented for constant factors!");

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
   static_assert(domain_of_integration != DomainOfIntegration::face
                 || type == IntegrationType::normalVector
                 || type == IntegrationType::normalSign,
                 "IntegrationType not implemented on boundary.");


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
    constexpr auto lhsType = (type == IntegrationType::valueValue ||
                              type == IntegrationType::valueGrad)
                             ? EvaluationType::value : EvaluationType::grad;

    std::vector<FieldVector<double,1> > lhsValues =
        detail::LocalFunctionEvaluation<dim, lhsType, domain_of_integration>()
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
        detail::LocalFunctionEvaluation<dim, rhsType, domain_of_integration>()
                      (rhsLocalFiniteElement,
                       quadPos,
                       geometry,
                       rhsBeta);

    // Compute the actual matrix entries
    for (size_t i=0; i<nLhs; i++)
    {
      for (size_t j=0; j<nRhs; j++)
      {
        elementMatrix[i+lhsSpaceOffset][j+rhsSpaceOffset]
                += (lhsValues[i] * rhsValues[j]) * factor
                   * quad[pt].weight() * integrationElement;
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
    const double integrationElement = intersection.geometry().integrationElement(quadFacePos);
    //const double integrationElement = integrationOuterNormal.two_norm();

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
        if(type == IntegrationType::normalVector) {
        elementMatrix[i+lhsSpaceOffset][j+rhsSpaceOffset]
                += ((lhsBeta*integrationOuterNormal) * factor
                    * lhsValues[i] * rhsValues[j]) * quadFace[pt].weight();
        } else if(type == IntegrationType::normalSign) {
        elementMatrix[i+lhsSpaceOffset][j+rhsSpaceOffset]
                += (sign * factor * lhsValues[i] * rhsValues[j])
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
