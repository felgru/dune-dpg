// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_INTEGRALTERM_HH
#define DUNE_DPG_INTEGRALTERM_HH

#include <tuple>

#include <dune/geometry/quadraturerules/splitquadraturerule.hh>
#include <dune/istl/matrix.hh>

#include "assemble_types.hh"
#include "type_traits.hh"
#include "quadrature.hh"
#include "localevaluation.hh"

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


namespace detail {
  template <IntegrationType type,
            class LhsSpace,
            class RhsSpace,
            bool = is_RefinedFiniteElement<LhsSpace>::value,
            bool = is_RefinedFiniteElement<RhsSpace>::value>
  struct GetLocalMatrix
  {
    using LhsLocalView = typename LhsSpace::LocalView;
    using RhsLocalView = typename RhsSpace::LocalView;

    template <class MatrixType,
              class Element,
              class FactorType,
              class DirectionType>
    inline static void interiorImpl(const LhsLocalView&,
                                    const RhsLocalView&,
                                    MatrixType&,
                                    size_t,
                                    size_t,
                                    unsigned int,
                                    const Element&,
                                    const FactorType&,
                                    const DirectionType&,
                                    const DirectionType&);

    template <class MatrixType,
              class Intersection,
              class FactorType,
              class DirectionType>
    inline static void faceImpl(const LhsLocalView&,
                                const RhsLocalView&,
                                MatrixType&,
                                size_t,
                                size_t,
                                unsigned int,
                                const Intersection&,
                                const FactorType&,
                                const DirectionType&,
                                const DirectionType&);
  };
}



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
         typename std::enable_if<
                     integrationType == IntegrationType::valueValue
                  || integrationType == IntegrationType::normalSign>::type*
                = nullptr
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
         typename std::enable_if<
                     integrationType == IntegrationType::gradValue
                  || integrationType == IntegrationType::valueGrad
                  || integrationType == IntegrationType::gradGrad
                  || integrationType == IntegrationType::normalVector
                  || integrationType == IntegrationType::travelDistanceWeighted
                  >::type*
           = nullptr
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
         typename std::enable_if<
                     integrationType == IntegrationType::gradGrad>::type*
                = nullptr
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
  static_assert(std::is_same<typename std::decay<DirectionType>::type,
                             FieldVector<double, 2>
                            >::value,
             "getLocalMatrix only implemented for constant flow!");

  static_assert(type == IntegrationType::valueValue
             || type == IntegrationType::gradValue
             || type == IntegrationType::valueGrad
             || type == IntegrationType::gradGrad
             || type == IntegrationType::normalVector
             || type == IntegrationType::normalSign
             || type == IntegrationType::travelDistanceWeighted,
             "Use of unknown IntegrationType.");
  static_assert(domain_of_integration != DomainOfIntegration::interior
                || type == IntegrationType::valueValue
                || type == IntegrationType::gradValue
                || type == IntegrationType::valueGrad
                || type == IntegrationType::gradGrad,
                "IntegrationType not implemented on interior.");
  static_assert(domain_of_integration != DomainOfIntegration::face
                || type == IntegrationType::normalVector
                || type == IntegrationType::normalSign
                || type == IntegrationType::travelDistanceWeighted,
                "IntegrationType not implemented on boundary.");

  using LhsSpace = typename LhsLocalView::GlobalBasis;
  using RhsSpace = typename RhsLocalView::GlobalBasis;

  // Get the grid element from the local FE basis view
  using Element = typename LhsLocalView::Element;
  const Element& element = lhsLocalView.element();

  const auto lhsOrder = lhsLocalView.tree().finiteElement().localBasis().order();
  const auto rhsOrder = rhsLocalView.tree().finiteElement().localBasis().order();

  /* TODO: We might need a higher order when factor is a function. */
  /* TODO: Assuming β const. */
  const unsigned int quadratureOrder = lhsOrder + rhsOrder;


  if(domain_of_integration == DomainOfIntegration::interior) {
    detail::GetLocalMatrix<type, LhsSpace, RhsSpace>
                         ::interiorImpl(lhsLocalView,
                                        rhsLocalView,
                                        elementMatrix,
                                        lhsSpaceOffset,
                                        rhsSpaceOffset,
                                        quadratureOrder,
                                        element,
                                        factor,
                                        lhsBeta,
                                        rhsBeta);
  } else {
    detail::GetLocalMatrix<type, LhsSpace, RhsSpace>
                         ::faceImpl(lhsLocalView,
                                    rhsLocalView,
                                    elementMatrix,
                                    lhsSpaceOffset,
                                    rhsSpaceOffset,
                                    quadratureOrder,
                                    element,
                                    factor,
                                    lhsBeta,
                                    rhsBeta);

  }
}

namespace detail {
  // This function expects data transformed to the reference triangle:
  // A point in the inflow boundary of the reference cell and the
  // transport direction transformed under the global-to-local mapping.
  template<class ReferenceCellCoordinate, class ReferenceCellDirection>
  double travelDistance(
      const ReferenceCellCoordinate& x,
      const ReferenceCellDirection& beta)
  {
    if(x[0]) { /* x[0] != 0 */
      if(x[1]) { /* x[1] != 0 */
        double a = x[1] - x[0]*beta[1]/beta[0];
        if(0 <= a && a <= 1)
          return -x[0]/beta[0];
        else
          return -x[1]/beta[1];
      } else { /* x[1] == 0 */
        double a = -beta[1]/beta[0]*x[0];
        if(0 <= a && a <= 1)
          return -x[0]/beta[0];
        else
          return (1-x[0])/(beta[0]+beta[1]);
      }
    } else { /* x[0] == 0 */
      double b = -beta[0]/beta[1]*x[1];
      if(0 <= b && b <= 1)
        return -x[1]/beta[1];
      else
        return (1-x[1])/(beta[0]+beta[1]);
    }
  }

  template<class Intersection, class ReferenceCellDirection>
  double splitPointOfInflowFaceInTriangle(
      const Intersection& intersection,
      ReferenceCellDirection& referenceBeta)
  {
    // This gets a bit ugly as we have to check the orientation of the face
    double splitPoint;
    const double tol = 1e-10;
    if(fabs(referenceBeta[0]) < tol) referenceBeta[0] = 0.;
    if(fabs(referenceBeta[1]) < tol) referenceBeta[1] = 0.;

    if(referenceBeta[0] > 0) {
      if((intersection.geometryInInside().global({0})
          - FieldVector<double,2>{0.,0.}).two_norm() < tol)
        splitPoint = -referenceBeta[1]/referenceBeta[0];
      else
        splitPoint = 1.+referenceBeta[1]/referenceBeta[0];
    } else if(referenceBeta[1] > 0) {
      if((intersection.geometryInInside().global({0})
          - FieldVector<double,2>{0.,0.}).two_norm() < tol)
        splitPoint = -referenceBeta[0]/referenceBeta[1];
      else
        splitPoint = 1.+referenceBeta[0]/referenceBeta[1];
    } else {
      if((intersection.geometryInInside().global({0})
          - FieldVector<double,2>{0.,1.}).two_norm() < tol)
        splitPoint = referenceBeta[0]/(referenceBeta[0]+referenceBeta[1]);
      else
        splitPoint = 1.-referenceBeta[0]/(referenceBeta[0]+referenceBeta[1]);
    }
    if(fabs(splitPoint)< tol)         splitPoint = 0.;
    else if(fabs(splitPoint-1.)< tol) splitPoint = 1.;

    assert(splitPoint >= 0 && splitPoint <= 1);
    return splitPoint;
  }

  template <class Geometry, class SubGeometry, int dim>
  FieldVector<double,dim> referenceBeta(
      const Geometry& geometry,
      const SubGeometry& subGeometryInReferenceElement,
      const FieldVector<double, dim>& beta)
  {
    static_assert(dim==2, "Computation of transport direction on reference"
                          " cell only implemented in 2d!");
    /* This won't work for curvilinear elements, but they don't seem
     * to be supported by UG anyway. */
    const auto& jacobianSubInverse
        = subGeometryInReferenceElement.jacobianInverseTransposed({0., 0.});
    const auto& jacobianInverse = geometry.jacobianInverseTransposed({0., 0.});
    FieldVector<double,dim> referenceBetaSub, referenceBeta;
    jacobianInverse.mtv(beta, referenceBeta);
    jacobianSubInverse.mtv(referenceBeta, referenceBetaSub);

    return referenceBetaSub;
  }
}


} // end namespace Dune

#include "integralterm_uu_impl.hh"
#include "integralterm_rr_impl.hh"
#include "integralterm_ru_impl.hh"

#endif // DUNE_DPG_INTEGRALTERM_HH
