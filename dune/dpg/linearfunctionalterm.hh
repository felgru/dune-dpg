// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_LINEARFUNCTIONALTERM_HH
#define DUNE_DPG_LINEARFUNCTIONALTERM_HH

#include <tuple>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/localevaluation.hh>
#include <dune/dpg/quadrature.hh>

#include <dune/istl/bvector.hh>

#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/functional/generation/make_fused_procedure.hpp>


namespace Dune {

  template <DomainOfIntegration domain_of_integration,
            class SolutionSpace,
            class FunctionalVector>
  class LinearFunctionalTerm
  {
    using SolutionLocalView = typename SolutionSpace::LocalView;
    using SolutionLocalIndexSet = typename SolutionSpace::LocalIndexSet;
  public:

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
        , solutionLocalIndexSet(solutionSpace.localIndexSet())
    {};

    /**
     * \brief Compute the rhs vector for a single element.
     *
     * The local integrals will be added with the given offsets
     * to \p elementVector.
     *
     * \note Not thread-safe!
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
    FunctionalVector functionalVector;
    mutable SolutionLocalView solutionLocalView;
    mutable SolutionLocalIndexSet solutionLocalIndexSet;
  };

/**
 * \brief Creates a tuple of a LinearFunctionalTerm and the index
 *        of the space involved.
 *
 * \param functionalVector  the coefficients describing the functional in
 *                          the given solutionSpace basis
 * \param solutionSpace  the space in which the functional is defined
 * \tparam spaceIndex the index of the test space
 * \tparam domainOfIntegration
 */
template<size_t spaceIndex,
         DomainOfIntegration domainOfIntegration,
         class SolutionSpace,
         class FunctionalVector>
auto make_LinearFunctionalTerm(const FunctionalVector& functionalVector,
                               const SolutionSpace& solutionSpace)
    -> std::tuple<std::integral_constant<size_t, spaceIndex>,
                  LinearFunctionalTerm<domainOfIntegration,
                                       SolutionSpace,
                                       FunctionalVector> >
{
  return std::make_tuple(
              std::integral_constant<size_t, spaceIndex>(),
              LinearFunctionalTerm<domainOfIntegration,
                                   SolutionSpace,
                                   FunctionalVector>(functionalVector,
                                                     solutionSpace));
}


namespace detail {
  template <class TestSpace,
            class SolutionSpace,
            bool = is_RefinedFiniteElement<TestSpace>::value,
            bool = is_RefinedFiniteElement<SolutionSpace>::value>
  struct ApplyLocalFunctional;
}

#include "applylocalfunctional_uu_impl.hh"
#include "applylocalfunctional_ru_impl.hh"


template <DomainOfIntegration domain_of_integration,
          class SolutionSpace,
          class FunctionalVector>
template <class LocalView,
          class VectorType>
void LinearFunctionalTerm<domain_of_integration,
                          SolutionSpace, FunctionalVector>::
getLocalVector(const LocalView& localView,
               VectorType& elementVector,
               size_t spaceOffset) const
{
  solutionLocalView.bind(localView.element());
  solutionLocalIndexSet.bind(solutionLocalView);

  // Now get the local contribution to the right-hand side vector
  static_assert(domain_of_integration == DomainOfIntegration::interior,
      "LinearFunctionalTerm only implemented for interior integrals.");
  detail::ApplyLocalFunctional
    < typename LocalView::GlobalBasis
    , SolutionSpace
    >::interiorImpl(
          localView,
          solutionLocalView,
          elementVector,
          spaceOffset,
          solutionLocalIndexSet,
          localView.element(),
          functionalVector);
}

} // end namespace Dune

#endif // DUNE_DPG_LINEARFUNCTIONALTERM_HH
