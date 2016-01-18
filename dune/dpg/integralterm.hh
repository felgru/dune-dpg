// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_INTEGRALTERM_HH
#define DUNE_DPG_INTEGRALTERM_HH

#include <tuple>
#include <vector>

#include <dune/istl/matrix.hh>

#include <dune/functions/functionspacebases/interpolate.hh>

#include "assemble_types.hh"
#include "type_traits.hh"
#include "quadrature.hh"

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
                  || integrationType == IntegrationType::normalVector>::type*
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

/* We need to make this a class, as partial specializations of
 * function templates are not allowed. */
template<int dim, EvaluationType type,
         DomainOfIntegration domain_of_integration,
         bool isDGRefined>
struct LocalRefinedFunctionEvaluation {

  template <class LocalFiniteElement, class Geometry>
  std::vector<FieldVector<double,1> >
  operator() (const LocalFiniteElement& localFiniteElement,
              unsigned int subElement,
              const FieldVector<double, dim>& quadPos,
              const Geometry& geometry,
              const Geometry& subGeometryInReferenceElement,
              const FieldVector<double, dim>& beta) const;
};

template<int dim, DomainOfIntegration domain_of_integration>
struct LocalRefinedFunctionEvaluation<dim, EvaluationType::value,
                               domain_of_integration, false> {

  template <class LocalFiniteElement, class Geometry>
  std::vector<FieldVector<double,1> > operator()
                      (const LocalFiniteElement& localFiniteElement,
                       unsigned int,
                       const FieldVector<double, dim>& quadPos,
                       const Geometry& geometry,
                       const Geometry& subGeometryInReferenceElement,
                       const FieldVector<double, dim>&) const
  {
    // values of the shape functions
    std::vector<FieldVector<double,1> > values;
    localFiniteElement.localBasis().evaluateFunction(quadPos, values);
    return values;
  }
};

template<int dim, DomainOfIntegration domain_of_integration>
struct LocalRefinedFunctionEvaluation<dim, EvaluationType::grad,
                               domain_of_integration, false> {

  template <class LocalFiniteElement, class Geometry>
  std::vector<FieldVector<double,1> > operator()
                      (const LocalFiniteElement& localFiniteElement,
                       unsigned int,
                       const FieldVector<double, dim> & quadPos,
                       const Geometry& geometry,
                       const Geometry& subGeometryInReferenceElement,
                       const FieldVector<double, dim>& beta) const
  {
    const auto& jacobianSub
        = subGeometryInReferenceElement.jacobianInverseTransposed(quadPos);
    const auto& jacobian = geometry.jacobianInverseTransposed
                           (subGeometryInReferenceElement.global(quadPos));
    // The gradients of the shape functions on the reference element
    std::vector<FieldMatrix<double,1,dim> > referenceGradients;
    localFiniteElement.localBasis()
            .evaluateJacobian(quadPos, referenceGradients);

    // Compute the shape function gradients on the real element
    std::vector<FieldVector<double, 1> >
            derivatives(referenceGradients.size());
    for (size_t i=0, i_max=referenceGradients.size(); i<i_max; i++)
    {
      FieldVector<double,dim> gradientRef, gradient;
      jacobianSub.mv(referenceGradients[i][0], gradientRef);
      jacobian.mv(gradientRef, gradient);
      derivatives[i] = beta * gradient;
    }

    return derivatives;
  }
};

template<int dim, DomainOfIntegration domain_of_integration>
struct LocalRefinedFunctionEvaluation<dim, EvaluationType::value,
                               domain_of_integration, true> {

  template <class LocalFiniteElement, class Geometry>
  std::vector<FieldVector<double,1> > operator()
                      (const LocalFiniteElement& localFiniteElement,
                       unsigned int subElement,
                       const FieldVector<double, dim>& quadPos,
                       const Geometry& geometry,
                       const Geometry& subGeometryInReferenceElement,
                       const FieldVector<double, dim>&) const
  {
    // values of the shape functions
    std::vector<FieldVector<double,1> > values;
    localFiniteElement.localBasis().evaluateFunction(subElement, quadPos, values);
    return values;
  }
};

template<int dim, DomainOfIntegration domain_of_integration>
struct LocalRefinedFunctionEvaluation<dim, EvaluationType::grad,
                               domain_of_integration, true> {

  template <class LocalFiniteElement, class Geometry>
  std::vector<FieldVector<double,1> > operator()
                      (const LocalFiniteElement& localFiniteElement,
                       unsigned int subElement,
                       const FieldVector<double, dim> & quadPos,
                       const Geometry& geometry,
                       const Geometry& subGeometryInReferenceElement,
                       const FieldVector<double, dim>& beta) const
  {
    const auto& jacobianSub
        = subGeometryInReferenceElement.jacobianInverseTransposed(quadPos);
    const auto& jacobian = geometry.jacobianInverseTransposed
                           (subGeometryInReferenceElement.global(quadPos));
    // The gradients of the shape functions on the reference element
    std::vector<FieldMatrix<double,1,dim> > referenceGradients;
    localFiniteElement.localBasis()
            .evaluateJacobian(subElement, quadPos, referenceGradients);

    // Compute the shape function gradients on the real element
    std::vector<FieldVector<double, 1> >
            derivatives(referenceGradients.size());
    for (size_t i=0, i_max=referenceGradients.size(); i<i_max; i++)
    {
      FieldVector<double,dim> gradientRef, gradient;
      jacobianSub.mv(referenceGradients[i][0], gradientRef);
      jacobian.mv(gradientRef, gradient);
      derivatives[i] = beta * gradient;
    }

    return derivatives;
  }
};

template<class FactorType, class PositionType,
         typename std::enable_if<std::is_arithmetic<FactorType>::value>
                              ::type* = nullptr >
inline double evaluateFactor(FactorType factor, PositionType)
{
  return factor;
}

template<class FactorType, class PositionType,
         typename std::enable_if<std::is_function<FactorType>::value>
                              ::type* = nullptr >
inline double evaluateFactor(FactorType factor, PositionType x)
{
  return factor(x);
}

template<class FactorType, class PositionType,
         typename std::enable_if<is_vector<FactorType>::value>
                              ::type* = nullptr >
inline double evaluateFactor(const FactorType& factor, PositionType x)
{
  return factor[x];
}

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

  using LhsSpace = typename LhsLocalView::GlobalBasis;
  using RhsSpace = typename RhsLocalView::GlobalBasis;

  // Get the grid element from the local FE basis view
  using Element = typename std::remove_pointer<LhsLocalView>::type::Element;
  const Element& element = lhsLocalView.element();

  const auto lhsOrder = lhsLocalView.tree().finiteElement().localBasis().order();
  const auto rhsOrder = rhsLocalView.tree().finiteElement().localBasis().order();

  /* TODO: We might need a higher order when factor is a function. */
  /* TODO: Assuming β const. */
  unsigned int quadratureOrder = lhsOrder + rhsOrder;


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


// UGGrid needs to be forward-declared in case that we are not using
// refined elements who would already have included the UG headers.
template<int dim> class UGGrid;
namespace detail {

#include "integralterm_uu_impl.hh"
#include "integralterm_rr_impl.hh"
#include "integralterm_ru_impl.hh"

} // end namespace detail

} // end namespace Dune

#endif // DUNE_DPG_INTEGRALTERM_HH
