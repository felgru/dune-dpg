// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_LINEARINTEGRALTERM_HH
#define DUNE_DPG_LINEARINTEGRALTERM_HH

#include <tuple>
#include <vector>

#include <dune/istl/matrix.hh>

#include "assemble_types.hh"
#include "type_traits.hh"
#include "quadrature.hh"
#include "localevaluation.hh"
#include "getvolumeterm_impl.hh"

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
  template <// IntegrationType type, TODO: reintroduce
            // DomainOfIntegration domain_of_integration, TODO: reintroduce?
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
     * \pre The localViews have to be bound to the same element.
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
 * \brief Creates a Tuple of an LinearIntegralTerm and the indices
 *        of both spaces involved.
 *
 * \param c  the factor with which we multiply the integrand
 * \tparam spaceIndex the index of the test space
 * \tparam FactorType  the type of the factor \p c
 */
template<size_t spaceIndex,
         // IntegrationType integrationType,
         // DomainOfIntegration domainOfIntegration,
         class FactorType>
auto make_LinearIntegralTerm(FactorType c)
    -> std::tuple<std::integral_constant<size_t, spaceIndex>,
                  LinearIntegralTerm<FactorType> >
{
  return std::tuple<std::integral_constant<size_t, spaceIndex>,
                LinearIntegralTerm<FactorType> >
         ({}, LinearIntegralTerm<FactorType>(c));
}

/**
 * \brief Creates a Tuple of an LinearIntegralTerm and the indices
 *        of both spaces involved.
 *
 * \param c     the factor with which we multiply the integrand
 * \param beta  the transport direction
 * \tparam spaceIndex     the index of the test space
 * \tparam FactorType     the type of the factor \p c
 * \tparam DirectionType  the type of the transport direction \p beta
 */
template<size_t spaceIndex,
         // IntegrationType integrationType,
         // DomainOfIntegration domainOfIntegration,
         class FactorType, class DirectionType>
auto make_LinearIntegralTerm(FactorType c, DirectionType beta)
    -> std::tuple<std::integral_constant<size_t, spaceIndex>,
                  LinearIntegralTerm<FactorType, DirectionType> >
{
  return std::tuple<std::integral_constant<size_t, spaceIndex>,
                LinearIntegralTerm<FactorType, DirectionType> >
         ({}, LinearIntegralTerm<FactorType, DirectionType>(c, beta));
}


template<class FactorType, class DirectionType>
template <class LocalView,
          class VectorType>
void LinearIntegralTerm<FactorType, DirectionType>
     ::getLocalVector(
        const LocalView& localView,
        VectorType& elementVector,
        size_t spaceOffset) const
{
  static_assert(std::is_same<typename std::decay<DirectionType>::type,
                             FieldVector<double, 2>
                            >::value,
             "getLocalVector only implemented for constant flow!");

  detail::GetVolumeTerm_Impl<LocalView, FactorType>
             ::getVolumeTerm(localView,
                             elementVector,
                             // SpaceOffset,
                             // quadratureOrder,
                             // element,
                             factor
                             //, beta
                            );

}


} // end namespace Dune

#endif // DUNE_DPG_LINEARINTEGRALTERM_HH
