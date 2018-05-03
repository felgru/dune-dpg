// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_INTEGRALTERM_HH
#define DUNE_DPG_INTEGRALTERM_HH

#include <tuple>

#include <dune/geometry/quadraturerules/splitquadraturerule.hh>
#include <dune/istl/matrix.hh>

#include "assemble_types.hh"
#include "localcoefficients.hh"
#include "localevaluation.hh"
#include "quadrature.hh"
#include "traveldistancenorm.hh"
#include "type_traits.hh"

namespace Dune {

  /**
   * \brief This class describes an integral term.
   *
   * This is the essential building block from which BilinearForm and
   * InnerProduct are built.
   *
   * \tparam integrationType  the form of the integrand, see #IntegrationType
   * \tparam domainOfIntegration  see #DomainOfIntegration
   * \tparam LocalCoefficients  a class that contains the factor with which
   *                    we multiply the integrand and the transport directions
   */
  template <IntegrationType type,
            DomainOfIntegration domain_of_integration,
            class LocalCoefficients>
  class IntegralTerm
  {
    using LocalFactor = typename LocalCoefficients::LocalFactor;
    using Element = typename LocalCoefficients::Element;
  public:

    IntegralTerm () = delete;

    /**
     * \brief constructor for IntegralTerm
     *
     * \note For your convenience, use make_IntegralTerm() instead.
     */
    IntegralTerm (LocalCoefficients&& localCoefficients)
        : localCoefficients_(localCoefficients)
    {};

    /**
     * \brief Compute the stiffness matrix for a single element.
     *
     * The local integrals will be added with the given offsets
     * to \p elementMatrix.
     *
     * \pre The localViews have to be bound to the same element as the
     *      IntegralTerm.
     *
     * \param[in]     lhsLocalView    local view of the left space
     * \param[in]     rhsLocalView    local view of the right space
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

    void bind(const Element& element)
    {
      localCoefficients_.bind(element);
    }

  private:
    LocalCoefficients localCoefficients_;
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
              class LocalCoefficients>
    inline static void interiorImpl(const LhsLocalView&,
                                    const RhsLocalView&,
                                    MatrixType&,
                                    size_t,
                                    size_t,
                                    unsigned int,
                                    const Element&,
                                    const LocalCoefficients&);

    template <class MatrixType,
              class Intersection,
              class LocalCoefficients>
    inline static void faceImpl(const LhsLocalView&,
                                const RhsLocalView&,
                                MatrixType&,
                                size_t,
                                size_t,
                                unsigned int,
                                const Intersection&,
                                const LocalCoefficients&);
  };
}



/**
 * \brief Creates a Tuple of an IntegralTerm and the indices
 *        of both spaces involved.
 *
 * \param c  the factor with which we multiply the integrand
 * \tparam lhsSpaceIndex the index of the left space
 * \tparam rhsSpaceIndex the index of the right space
 * \tparam integrationType  the form of the integrand, see #IntegrationType
 * \tparam domainOfIntegration  see #DomainOfIntegration
 * \tparam Factor  the type of the factor \p c
 */
template<size_t lhsSpaceIndex,
         size_t rhsSpaceIndex,
         IntegrationType integrationType,
         DomainOfIntegration domainOfIntegration,
         class Factor,
         typename std::enable_if<
                     integrationType == IntegrationType::valueValue
                  || integrationType == IntegrationType::normalSign>::type*
                = nullptr
        >
auto make_IntegralTerm(Factor c)
    -> std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                  std::integral_constant<size_t, rhsSpaceIndex>,
                  IntegralTerm<integrationType, domainOfIntegration,
                               detail::LocalCoefficients::OnlyFactor<Factor>> >
{
  return std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                std::integral_constant<size_t, rhsSpaceIndex>,
                IntegralTerm<integrationType, domainOfIntegration,
                             detail::LocalCoefficients::OnlyFactor<Factor>> >
         ({},{},
          {detail::LocalCoefficients::OnlyFactor<Factor>(c)});
}

/**
 * \brief Creates a Tuple of an IntegralTerm and the indices
 *        of both spaces involved.
 *
 * \param c     the factor with which we multiply the integrand
 * \param beta  the transport direction
 * \tparam lhsSpaceIndex the index of the left space
 * \tparam rhsSpaceIndex the index of the right space
 * \tparam integrationType  the form of the integrand, see #IntegrationType
 * \tparam domainOfIntegration  see #DomainOfIntegration
 * \tparam Factor     the type of the factor \p c
 * \tparam Direction  the type of the transport direction \p beta
 */
template<size_t lhsSpaceIndex,
         size_t rhsSpaceIndex,
         IntegrationType integrationType,
         DomainOfIntegration domainOfIntegration,
         class Factor, class Direction,
         typename std::enable_if<
                     integrationType == IntegrationType::gradValue
                  || integrationType == IntegrationType::valueGrad
                  || integrationType == IntegrationType::gradGrad
                  || integrationType == IntegrationType::normalVector
                  || integrationType == IntegrationType::travelDistanceWeighted
                  >::type*
           = nullptr
        >
auto make_IntegralTerm(Factor c, Direction beta)
    -> std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                  std::integral_constant<size_t, rhsSpaceIndex>,
                  IntegralTerm<integrationType, domainOfIntegration,
                               detail::LocalCoefficients::FactorAndDirection
                                   <Factor, Direction>> >
{
  return std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                std::integral_constant<size_t, rhsSpaceIndex>,
                IntegralTerm<integrationType, domainOfIntegration,
                             detail::LocalCoefficients::FactorAndDirection
                                 <Factor, Direction>> >
         ({},{},
          {detail::LocalCoefficients::FactorAndDirection
                                 <Factor, Direction>(c, beta)});
}

/**
 * \brief Creates a Tuple of an IntegralTerm and the indices
 *        of both spaces involved.
 *
 * \param c        the factor with which we multiply the integrand
 * \param lhsBeta  the transport direction for the left space
 * \param rhsBeta  the transport direction for the right space
 * \tparam lhsSpaceIndex the index of the left space
 * \tparam rhsSpaceIndex the index of the right space
 * \tparam integrationType  the form of the integrand, see #IntegrationType
 * \tparam domainOfIntegration  see #DomainOfIntegration
 * \tparam Factor     the type of the factor \p c
 * \tparam Direction  the type of the transport directions
 */
template<size_t lhsSpaceIndex,
         size_t rhsSpaceIndex,
         IntegrationType integrationType,
         DomainOfIntegration domainOfIntegration,
         class Factor, class Direction,
         typename std::enable_if<
                     integrationType == IntegrationType::gradGrad>::type*
                = nullptr
        >
auto make_IntegralTerm(Factor c,
                       Direction lhsBeta,
                       Direction rhsBeta)
    -> std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                  std::integral_constant<size_t, rhsSpaceIndex>,
                  IntegralTerm<integrationType, domainOfIntegration,
                               detail::LocalCoefficients::FactorAndTwoDirections
                                 <Factor, Direction, Direction>> >
{
  return std::tuple<std::integral_constant<size_t, lhsSpaceIndex>,
                std::integral_constant<size_t, rhsSpaceIndex>,
                IntegralTerm<integrationType, domainOfIntegration,
                             detail::LocalCoefficients::FactorAndTwoDirections
                                 <Factor, Direction, Direction>> >
         ({},{},
          {detail::LocalCoefficients::FactorAndTwoDirections
                   <Factor, Direction, Direction>(c, lhsBeta, rhsBeta)});
}


template<IntegrationType type, DomainOfIntegration domain_of_integration,
         class LocalCoefficients>
template <class LhsLocalView,
          class RhsLocalView,
          class MatrixType>
void IntegralTerm<type, domain_of_integration, LocalCoefficients>
     ::getLocalMatrix(
        const LhsLocalView& lhsLocalView,
        const RhsLocalView& rhsLocalView,
        MatrixType& elementMatrix,
        const size_t lhsSpaceOffset,
        const size_t rhsSpaceOffset) const
{
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
  const auto quadratureOrder = lhsOrder + rhsOrder;


  using namespace Dune::Hybrid;
  ifElse(equals(std::integral_constant<DomainOfIntegration,
                                       domain_of_integration>(),
                std::integral_constant<DomainOfIntegration,
                                       DomainOfIntegration::interior>()),
  [&](auto id) {
    detail::GetLocalMatrix<type, LhsSpace, RhsSpace>
                         ::interiorImpl(lhsLocalView,
                                        rhsLocalView,
                                        elementMatrix,
                                        lhsSpaceOffset,
                                        rhsSpaceOffset,
                                        quadratureOrder,
                                        element,
                                        localCoefficients_);
  }, [&](auto id) {
    detail::GetLocalMatrix<type, LhsSpace, RhsSpace>
                         ::faceImpl(lhsLocalView,
                                    rhsLocalView,
                                    elementMatrix,
                                    lhsSpaceOffset,
                                    rhsSpaceOffset,
                                    quadratureOrder,
                                    element,
                                    localCoefficients_);

  });
}

} // end namespace Dune

#include "integralterm_uu_impl.hh"
#include "integralterm_rr_impl.hh"
#include "integralterm_ru_impl.hh"

#endif // DUNE_DPG_INTEGRALTERM_HH
