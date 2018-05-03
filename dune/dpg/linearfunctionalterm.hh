// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_LINEARFUNCTIONALTERM_HH
#define DUNE_DPG_LINEARFUNCTIONALTERM_HH

#include <tuple>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/localevaluation.hh>
#include <dune/dpg/quadrature.hh>

#include <dune/istl/bvector.hh>


namespace Dune {

  /**
   * \brief This class describes an integration of a test function against
   *        a function given as an FE coefficient vector.
   *
   * This is one of the essential building blocks from which a LinearForm
   * is built. It implements
   * \f[ \int_\Omega v g \,\mathrm dx \f]
   * where \f$ g \f$ is a function given by a SolutionSpace and a vector of
   * coefficients. \f$ v \f$ are the test functions.
   *
   * \tparam SolutionSpace  the GlobalBasis for the function against which
   *                        the test functions get integrated
   * \tparam FunctionalVector  the coefficient vector for aforementioned
   *                           function
   */
  template <class SolutionSpace,
            class FunctionalVector>
  class LinearFunctionalTerm
  {
    using SolutionLocalView = typename SolutionSpace::LocalView;
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    using SolutionLocalIndexSet = typename SolutionSpace::LocalIndexSet;
#endif
  public:
    using Element = typename SolutionLocalView::Element;

    LinearFunctionalTerm () = delete;

    /**
     * \brief constructor for LinearFunctionalTerm
     *
     * \note For your convenience, use make_LinearFunctionalTerm() instead.
     */
    LinearFunctionalTerm (const FunctionalVector& functionalVector,
                          const SolutionSpace& solutionSpace)
        : functionalVector(functionalVector)
        , solutionLocalView(solutionSpace.localView())
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
        , solutionLocalIndexSet(solutionSpace.localIndexSet())
#endif
    {};

    /**
     * \brief Compute the rhs vector for a single element.
     *
     * The local integrals will be added with the given offsets
     * to \p elementVector.
     *
     * \pre The localView has to be bound to the same element as the
     *      LinearFunctionalTerm.
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

    void bind(const Element& element)
    {
      solutionLocalView.bind(element);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
      solutionLocalIndexSet.bind(solutionLocalView);
#endif
    }

  private:
    const FunctionalVector& functionalVector;
    SolutionLocalView solutionLocalView;
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    SolutionLocalIndexSet solutionLocalIndexSet;
#endif
  };

/**
 * \brief Creates a tuple of a LinearFunctionalTerm and the index
 *        of the space involved.
 *
 * \param functionalVector  the coefficients describing the functional in
 *                          the given solutionSpace basis
 * \param solutionSpace  the space in which the functional is defined
 * \tparam spaceIndex the index of the test space
 * \tparam domainOfIntegration  see #DomainOfIntegration,
 *                              only interior is allowed here
 */
template<size_t spaceIndex,
         DomainOfIntegration domainOfIntegration,
         class SolutionSpace,
         class FunctionalVector>
[[deprecated("make_LinearFunctionalTerm with DomainOfIntegration parameter "
             "has been deprecated.")]]
auto make_LinearFunctionalTerm(const FunctionalVector& functionalVector,
                               const SolutionSpace& solutionSpace)
    -> std::tuple<std::integral_constant<size_t, spaceIndex>,
                  LinearFunctionalTerm<SolutionSpace,
                                       FunctionalVector> >
{
  static_assert(domainOfIntegration == DomainOfIntegration::interior,
                "make_LinearFunctionalTerm only implemented for interior "
                "domain of integration.");
  return std::make_tuple(
              std::integral_constant<size_t, spaceIndex>(),
              LinearFunctionalTerm<SolutionSpace,
                                   FunctionalVector>(functionalVector,
                                                     solutionSpace));
}

/**
 * \brief Creates a tuple of a LinearFunctionalTerm and the index
 *        of the space involved.
 *
 * \param functionalVector  the coefficients describing the functional in
 *                          the given solutionSpace basis
 * \param solutionSpace  the space in which the functional is defined
 * \tparam spaceIndex the index of the test space
 */
template<size_t spaceIndex,
         class SolutionSpace,
         class FunctionalVector>
auto make_LinearFunctionalTerm(const FunctionalVector& functionalVector,
                               const SolutionSpace& solutionSpace)
    -> std::tuple<std::integral_constant<size_t, spaceIndex>,
                  LinearFunctionalTerm<SolutionSpace,
                                       FunctionalVector> >
{
  return std::make_tuple(
              std::integral_constant<size_t, spaceIndex>(),
              LinearFunctionalTerm<SolutionSpace,
                                   FunctionalVector>(functionalVector,
                                                     solutionSpace));
}


namespace detail {
  template <IntegrationType type,
            class TestSpace,
            class SolutionSpace,
            bool = is_RefinedFiniteElement<TestSpace>::value,
            bool = is_RefinedFiniteElement<SolutionSpace>::value>
  struct ApplyLocalFunctional;
}

} // End namespace Dune

#include "applylocalfunctional_uu_impl.hh"
#include "applylocalfunctional_ru_impl.hh"
#include "applylocalfunctional_rr_impl.hh"

namespace Dune {

template <class SolutionSpace,
          class FunctionalVector>
template <class LocalView,
          class VectorType>
void LinearFunctionalTerm<SolutionSpace, FunctionalVector>::
getLocalVector(const LocalView& localView,
               VectorType& elementVector,
               const size_t spaceOffset) const
{
  // Now get the local contribution to the right-hand side vector
  detail::ApplyLocalFunctional
    < IntegrationType::valueValue
    , typename LocalView::GlobalBasis
    , SolutionSpace
    >::interiorImpl(
          localView,
          solutionLocalView,
          elementVector,
          spaceOffset,
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
          solutionLocalIndexSet,
#endif
          localView.element(),
          functionalVector);
}

  /**
   * \brief This class describes an integration of a test function against
   *        a function given as an FE coefficient vector over the skeleton
   *        of the grid.
   *
   * This is one of the essential building blocks from which a LinearForm
   * is built.
   *
   * \tparam integrationType  the form of the integrand, see #IntegrationType
   * \tparam SolutionSpace  the GlobalBasis for the function against which
   *                        the test functions get integrated
   * \tparam FunctionalVector  the coefficient vector for aforementioned
   *                           function
   * \tparam Factor     the type of the factor with which
   *                        we multiply the integrand
   * \tparam Direction  the type of the transport directions
   */
  template <IntegrationType type,
            class SolutionSpace,
            class FunctionalVector,
            class Factor,
            class Direction>
  class SkeletalLinearFunctionalTerm
  {
    using SolutionLocalView = typename SolutionSpace::LocalView;
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    using SolutionLocalIndexSet = typename SolutionSpace::LocalIndexSet;
#endif
    using LocalFactor = typename Factor::LocalFunction;
  public:
    using Element = typename SolutionLocalView::Element;

    SkeletalLinearFunctionalTerm () = delete;

    /**
     * \brief constructor for SkeletalLinearFunctionalTerm
     *
     * \note For your convenience, use make_SkeletalLinearFunctionalTerm()
     *       instead.
     */
    SkeletalLinearFunctionalTerm (const FunctionalVector& functionalVector,
                                  const SolutionSpace& solutionSpace,
                                  Factor coefficient,
                                  Direction beta)
        : factor(coefficient)
        , localFactor(localFunction(factor))
        , beta(beta)
        , functionalVector(functionalVector)
        , solutionLocalView(solutionSpace.localView())
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
        , solutionLocalIndexSet(solutionSpace.localIndexSet())
#endif
    {};

    /**
     * \brief Compute the rhs vector for a single element.
     *
     * The local integrals will be added with the given offsets
     * to \p elementVector.
     *
     * \pre The localView has to be bound to the same element as the
     *      LinearFunctionalTerm.
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

    void bind(const Element& element)
    {
      solutionLocalView.bind(element);
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
      solutionLocalIndexSet.bind(solutionLocalView);
#endif
      localFactor.bind(element);
    }

  private:
    Factor factor;
    LocalFactor localFactor;
    Direction beta;
    const FunctionalVector& functionalVector;
    SolutionLocalView solutionLocalView;
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
    SolutionLocalIndexSet solutionLocalIndexSet;
#endif
  };

/**
 * \brief Creates a tuple of a SkeletalLinearFunctionalTerm and the index
 *        of the space involved.
 *
 * \param functionalVector  the coefficients describing the functional in
 *                          the given solutionSpace basis
 * \param solutionSpace  the space in which the functional is defined
 * \param c     the factor with which we multiply the integrand
 * \param beta  the transport direction
 * \tparam spaceIndex the index of the test space
 * \tparam integrationType  the form of the integrand, see #IntegrationType
 */
template<size_t spaceIndex,
         IntegrationType integrationType,
         class SolutionSpace,
         class FunctionalVector,
         class Factor, class Direction,
         typename std::enable_if<
                     integrationType == IntegrationType::normalVector
                  || integrationType == IntegrationType::travelDistanceWeighted
                  >::type*
           = nullptr
        >
auto make_SkeletalLinearFunctionalTerm(
    const FunctionalVector& functionalVector,
    const SolutionSpace& solutionSpace,
    Factor c, Direction beta)
    -> std::tuple<std::integral_constant<size_t, spaceIndex>,
                  SkeletalLinearFunctionalTerm<integrationType,
                                               SolutionSpace,
                                               FunctionalVector,
                                               Factor,
                                               Direction> >
{
  return std::make_tuple(
              std::integral_constant<size_t, spaceIndex>(),
              SkeletalLinearFunctionalTerm<integrationType,
                                           SolutionSpace,
                                           FunctionalVector,
                                           Factor,
                                           Direction>
                (functionalVector, solutionSpace, c, beta));
}


template <IntegrationType type,
          class SolutionSpace,
          class FunctionalVector,
          class Factor,
          class Direction>
template <class LocalView,
          class VectorType>
void SkeletalLinearFunctionalTerm<type, SolutionSpace, FunctionalVector,
                                  Factor, Direction>::
getLocalVector(const LocalView& localView,
               VectorType& elementVector,
               size_t spaceOffset) const
{
  // Now get the local contribution to the right-hand side vector
  detail::ApplyLocalFunctional
    < type
    , typename LocalView::GlobalBasis
    , SolutionSpace
    >::faceImpl(
          localView,
          solutionLocalView,
          elementVector,
          spaceOffset,
#if not(DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7))
          solutionLocalIndexSet,
#endif
          localView.element(),
          functionalVector,
          localFactor,
          beta);
}

} // end namespace Dune

#endif // DUNE_DPG_LINEARFUNCTIONALTERM_HH
