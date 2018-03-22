// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_LINEARINTEGRALTERM_HH
#define DUNE_DPG_LINEARINTEGRALTERM_HH

#include <tuple>
#include <vector>

#include "assemble_types.hh"
#include "type_traits.hh"
#include "quadrature.hh"
#include "localevaluation.hh"
#include "locallinearterm_impl.hh"

namespace Dune {

  /**
   * \brief This class describes an integral term.
   *
   * This is the essential building block from which each LinearForm is built.
   *
   * \tparam integrationType  the form of the integrand, see #IntegrationType
   * \tparam domainOfIntegration  see #DomainOfIntegration
   * \tparam FactorType     the type of the factor with which
   *                        we multiply the integrand
   * \tparam DirectionType  the type of the transport directions
   */
  template <LinearIntegrationType type,
            DomainOfIntegration domain_of_integration,
            class FactorType,
            class DirectionType = FieldVector<double, 2> >
  class LinearIntegralTerm
  {
  public:

    LinearIntegralTerm () = delete;

    /**
     * \brief constructor for LinearIntegralTerm
     *
     * \note For your convenience, use make_LinearIntegralTerm() instead.
     */
    LinearIntegralTerm (FactorType factor = 1,
                        DirectionType beta = {1,1})
        : factor(factor),
          beta(beta)
    {};

    /**
     * \brief Compute the rhs vector for a single element.
     *
     * The local integrals will be added with the given offsets
     * to \p elementVector.
     *
     * \param[in]     localView       local view of the test space
     * \param[in,out] elementVector   the local rhs vector
     * \param         spaceOffset     row offset for the test space
     */
    template <class LocalView,
              class VectorType>
    void getLocalVector(const LocalView& localView,
                        VectorType& elementVector,
                        size_t spaceOffset) const;

  private:
    FactorType factor;
    DirectionType beta;

  };


/**
 * \brief Creates a tuple of a LinearIntegralTerm and the index
 *        of the space involved.
 *
 * \param c  the factor with which we multiply the integrand
 * \tparam spaceIndex the index of the test space
 * \tparam integrationType  the form of the integrand,
 *         see #LinearIntegrationType
 * \tparam domainOfIntegration  see #DomainOfIntegration
 * \tparam FactorType  the type of the factor \p c
 */
template<size_t spaceIndex,
         LinearIntegrationType integrationType,
         DomainOfIntegration domainOfIntegration,
         class FactorType>
auto make_LinearIntegralTerm(FactorType c)
    -> std::tuple<std::integral_constant<size_t, spaceIndex>,
                  LinearIntegralTerm<integrationType,
                                     domainOfIntegration,
                                     FactorType> >
{
  return std::make_tuple(
              std::integral_constant<size_t, spaceIndex>(),
              LinearIntegralTerm<integrationType,
                                 domainOfIntegration,
                                 FactorType>(c));
}

/**
 * \brief Creates a tuple of a LinearIntegralTerm and the index
 *        of the space involved.
 *
 * \param c     the factor with which we multiply the integrand
 * \param beta  the transport direction
 * \tparam spaceIndex     the index of the test space
 * \tparam integrationType  the form of the integrand,
 *         see #LinearIntegrationType
 * \tparam domainOfIntegration  see #DomainOfIntegration
 * \tparam FactorType     the type of the factor \p c
 * \tparam DirectionType  the type of the transport direction \p beta
 */
template<size_t spaceIndex,
         LinearIntegrationType integrationType,
         DomainOfIntegration domainOfIntegration,
         class FactorType,
         class DirectionType>
auto make_LinearIntegralTerm(FactorType c, DirectionType beta)
    -> std::tuple<std::integral_constant<size_t, spaceIndex>,
                  LinearIntegralTerm<integrationType,
                                     domainOfIntegration,
                                     FactorType,
                                     DirectionType> >
{
  return std::make_tuple(
              std::integral_constant<size_t, spaceIndex>(),
              LinearIntegralTerm<integrationType,
                                 domainOfIntegration,
                                 FactorType,
                                 DirectionType>(c, beta));
}


template<LinearIntegrationType integrationType,
         DomainOfIntegration domainOfIntegration,
         class FactorType,
         class DirectionType>
template <class LocalView,
          class VectorType>
void LinearIntegralTerm<integrationType,
                        domainOfIntegration,
                        FactorType,
                        DirectionType>
     ::getLocalVector(
        const LocalView& localView,
        VectorType& elementVector,
        const size_t spaceOffset) const
{
  static_assert(std::is_same<typename std::decay<DirectionType>::type,
                             FieldVector<double, 2>
                            >::value,
                "getLocalVector only implemented for constant flow!");

  static_assert(domainOfIntegration == DomainOfIntegration::interior,
                "DomainOfIntegration not implemented.");  //TODO

  static_assert(integrationType == LinearIntegrationType::valueFunction
             || integrationType == LinearIntegrationType::gradFunction,
                "LinearIntegrationType not implemented.");  //TODO

  using Space = typename LocalView::GlobalBasis;

  /* TODO: We might need a higher order when factor is a function. */
  /* TODO: Assuming β const. */
  const auto quadratureOrder = localView.tree().finiteElement().localBasis().order();
  if(domainOfIntegration == DomainOfIntegration::interior) {
    detail::GetLocalLinearTermVector<integrationType, Space>
               ::getLocalVector(localView,
                                elementVector,
                                spaceOffset,
                                quadratureOrder,
                                factor,
                                beta
                                );
  } else {
  //TODO
  }

}


} // end namespace Dune

#endif // DUNE_DPG_LINEARINTEGRALTERM_HH
