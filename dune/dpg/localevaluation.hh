// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_LOCALEVALUATION_HH
#define DUNE_DPG_LOCALEVALUATION_HH

#include <tuple>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

#include "assemble_types.hh"
#include "type_traits.hh"
#include <dune/common/std/type_traits.hh>

namespace Dune {
namespace detail {

/* We need to make this a class, as partial specializations of
 * function templates are not allowed. */
template<int dim, class LocalFiniteElement>
std::vector<FieldVector<double,1> >
evaluateLocalFunctionValue
                      (const LocalFiniteElement& localFiniteElement,
                       const FieldVector<double, dim>& quadPos)
{
  // values of the shape functions
  std::vector<FieldVector<double,1> > values;
  localFiniteElement.localBasis().evaluateFunction(quadPos, values);
  return values;
}

template<int dim, class LocalFiniteElement, class Geometry>
std::vector<FieldVector<double,1> >
evaluateLocalFunctionGrad
                      (const LocalFiniteElement& localFiniteElement,
                       const FieldVector<double, dim> & quadPos,
                       const Geometry& geometry,
                       const FieldVector<double, dim>& beta)
{
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
}

template<int dim, IntegrationType integrationType>
struct LocalFunctionEvaluation {
  constexpr static auto lhsType
    = (integrationType == IntegrationType::valueValue ||
       integrationType == IntegrationType::valueGrad)
      ? EvaluationType::value : EvaluationType::grad;

  constexpr static auto rhsType
    = (integrationType == IntegrationType::valueValue ||
       integrationType == IntegrationType::gradValue)
      ? EvaluationType::value : EvaluationType::grad;

  template<class LocalFiniteElement, class Geometry, class LocalCoefficients>
  static std::vector<FieldVector<double,1> >
  evaluateLhs(const LocalFiniteElement& localFiniteElement,
              const FieldVector<double, dim>& quadPos,
              const Geometry& geometry,
              const LocalCoefficients& localCoefficients)
  {
    return evaluateLhs_<lhsType>
                       (localFiniteElement,
                        quadPos,
                        geometry,
                        localCoefficients);
  }

  template<class LocalFiniteElement, class Geometry, class LocalCoefficients>
  static std::vector<FieldVector<double,1> >
  evaluateRhs(const LocalFiniteElement& localFiniteElement,
              const FieldVector<double, dim>& quadPos,
              const Geometry& geometry,
              const LocalCoefficients& localCoefficients)
  {
    return evaluateRhs_<rhsType>
                       (localFiniteElement,
                        quadPos,
                        geometry,
                        localCoefficients);
  }

  private:

  template<EvaluationType type,
           class LocalFiniteElement, class Geometry, class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::value>* = nullptr>
  static std::vector<FieldVector<double,1> >
  evaluateLhs_(const LocalFiniteElement& localFiniteElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& ,
               const LocalCoefficients& ) {
    return evaluateLocalFunctionValue(localFiniteElement, quadPos);
  }

  template<EvaluationType type,
           class LocalFiniteElement, class Geometry, class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::grad>* = nullptr>
  static std::vector<FieldVector<double,1> >
  evaluateLhs_(const LocalFiniteElement& localFiniteElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& geometry,
               const LocalCoefficients& localCoefficients) {
    return evaluateLocalFunctionGrad
                    (localFiniteElement,
                     quadPos,
                     geometry,
                     localCoefficients.localDirection()(quadPos));
  }

  template<EvaluationType type,
           class LocalFiniteElement, class Geometry, class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::value>* = nullptr>
  static std::vector<FieldVector<double,1> >
  evaluateRhs_(const LocalFiniteElement& localFiniteElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& ,
               const LocalCoefficients& ) {
    return evaluateLocalFunctionValue(localFiniteElement, quadPos);
  }

  template<EvaluationType type,
           class LocalFiniteElement, class Geometry, class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::grad>* = nullptr>
  static std::vector<FieldVector<double,1> >
  evaluateRhs_(const LocalFiniteElement& localFiniteElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& geometry,
               const LocalCoefficients& localCoefficients) {
    return evaluateLocalFunctionGrad
                    (localFiniteElement,
                     quadPos,
                     geometry,
                     localCoefficients.localSecondDirection()(quadPos));
  }
};

template<int dim, LinearIntegrationType integrationType>
struct LocalLinearTermFunctionEvaluation {
  constexpr static auto evaluationType = (integrationType ==
                      LinearIntegrationType::valueFunction)
                      ? EvaluationType::value : EvaluationType::grad;

  template<class LocalFiniteElement, class Geometry, class LocalCoefficients,
           std::enable_if_t<evaluationType
                            == EvaluationType::value>* = nullptr>
  static std::vector<FieldVector<double,1> >
  evaluate(const LocalFiniteElement& localFiniteElement,
           const FieldVector<double, dim>& quadPos,
           const Geometry& ,
           const LocalCoefficients& ) {
    return evaluateLocalFunctionValue(localFiniteElement, quadPos);
  }

  template<class LocalFiniteElement, class Geometry, class LocalCoefficients,
           std::enable_if_t<evaluationType
                            == EvaluationType::grad>* = nullptr>
  static std::vector<FieldVector<double,1> >
  evaluate(const LocalFiniteElement& localFiniteElement,
           const FieldVector<double, dim>& quadPos,
           const Geometry& geometry,
           const LocalCoefficients& localCoefficients) {
    return evaluateLocalFunctionGrad
                    (localFiniteElement,
                     quadPos,
                     geometry,
                     localCoefficients.localDirection()(quadPos));
  }
};

/* We need to make this a class, as partial specializations of
 * function templates are not allowed. */
template<bool isDGRefined>
struct LocalRefinedFunctionEvaluationHelper {

  template <int dim, class LocalFiniteElement,
            class Geometry, class SubGeometry>
  static std::vector<FieldVector<double,1> >
  evaluateValue(const LocalFiniteElement& localFiniteElement,
                unsigned int subElement,
                const FieldVector<double, dim>& quadPos,
                const Geometry& geometry,
                const SubGeometry& subGeometryInReferenceElement,
                const FieldVector<double, dim>& beta);

  template <int dim, class LocalFiniteElement,
            class Geometry, class SubGeometry>
  static std::vector<FieldVector<double,1> >
  evaluateGrad(const LocalFiniteElement& localFiniteElement,
               unsigned int subElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& geometry,
               const SubGeometry& subGeometryInReferenceElement,
               const FieldVector<double, dim>& beta);
};

template<>
struct LocalRefinedFunctionEvaluationHelper<false> {

  template<int dim, class LocalFiniteElement>
  static std::vector<FieldVector<double,1> >
  evaluateValue(const LocalFiniteElement& localFiniteElement,
                unsigned int,
                const FieldVector<double, dim>& quadPos)
  {
    // values of the shape functions
    std::vector<FieldVector<double,1> > values;
    localFiniteElement.localBasis().evaluateFunction(quadPos, values);
    return values;
  }

  template<int dim, class LocalFiniteElement,
           class Geometry, class SubGeometry>
  static std::vector<FieldVector<double,1> >
  evaluateGrad(const LocalFiniteElement& localFiniteElement,
               unsigned int,
               const FieldVector<double, dim> & quadPos,
               const Geometry& geometry,
               const SubGeometry& subGeometryInReferenceElement,
               const FieldVector<double, dim>& beta)
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

template<>
struct LocalRefinedFunctionEvaluationHelper<true> {

  template<int dim, class LocalFiniteElement>
  static std::vector<FieldVector<double,1> >
  evaluateValue(const LocalFiniteElement& localFiniteElement,
                unsigned int subElement,
                const FieldVector<double, dim>& quadPos)
  {
    // values of the shape functions
    std::vector<FieldVector<double,1> > values;
    localFiniteElement.localBasis().evaluateFunction(subElement, quadPos,
                                                     values);
    return values;
  }

  template<int dim, class LocalFiniteElement,
           class Geometry, class SubGeometry>
  static std::vector<FieldVector<double,1> >
  evaluateGrad(const LocalFiniteElement& localFiniteElement,
               unsigned int subElement,
               const FieldVector<double, dim> & quadPos,
               const Geometry& geometry,
               const SubGeometry& subGeometryInReferenceElement,
               const FieldVector<double, dim>& beta)
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

template<int dim, IntegrationType integrationType>
struct LocalRefinedFunctionEvaluation {
  constexpr static auto lhsType
    = (integrationType == IntegrationType::valueValue ||
       integrationType == IntegrationType::valueGrad)
      ? EvaluationType::value : EvaluationType::grad;

  constexpr static auto rhsType
    = (integrationType == IntegrationType::valueValue ||
       integrationType == IntegrationType::gradValue)
      ? EvaluationType::value : EvaluationType::grad;

  template<bool isDGRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients>
  static std::vector<FieldVector<double,1> >
  evaluateLhs(const LocalFiniteElement& localFiniteElement,
              unsigned int subElement,
              const FieldVector<double, dim>& quadPos,
              const Geometry& geometry,
              const SubGeometry& subGeometryInReferenceElement,
              const LocalCoefficients& localCoefficients)
  {
    return evaluateLhs_<lhsType, isDGRefined>
            (localFiniteElement,
             subElement,
             quadPos,
             geometry,
             subGeometryInReferenceElement,
             localCoefficients);
  }

  template<bool isDGRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients>
  static std::vector<FieldVector<double,1> >
  evaluateRhs(const LocalFiniteElement& localFiniteElement,
              unsigned int subElement,
              const FieldVector<double, dim>& quadPos,
              const Geometry& geometry,
              const SubGeometry& subGeometryInReferenceElement,
              const LocalCoefficients& localCoefficients)
  {
    return evaluateRhs_<rhsType, isDGRefined>
            (localFiniteElement,
             subElement,
             quadPos,
             geometry,
             subGeometryInReferenceElement,
             localCoefficients);
  }

  private:

  template<EvaluationType type, bool isDGRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::value>* = nullptr>
  static std::vector<FieldVector<double,1> >
  evaluateLhs_(const LocalFiniteElement& localFiniteElement,
               unsigned int subElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& ,
               const SubGeometry& ,
               const LocalCoefficients& ) {
    return LocalRefinedFunctionEvaluationHelper<isDGRefined>::
              evaluateValue(localFiniteElement, subElement, quadPos);
  }

  template<EvaluationType type, bool isDGRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::grad>* = nullptr>
  static std::vector<FieldVector<double,1> >
  evaluateLhs_(const LocalFiniteElement& localFiniteElement,
               unsigned int subElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& geometry,
               const SubGeometry& subGeometryInReferenceElement,
               const LocalCoefficients& localCoefficients) {
    return LocalRefinedFunctionEvaluationHelper<isDGRefined>::
              evaluateGrad(localFiniteElement,
                           subElement,
                           quadPos,
                           geometry,
                           subGeometryInReferenceElement,
                           localCoefficients.localDirection()(quadPos));
  }

  template<EvaluationType type, bool isDGRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::value>* = nullptr>
  static std::vector<FieldVector<double,1> >
  evaluateRhs_(const LocalFiniteElement& localFiniteElement,
               unsigned int subElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& ,
               const SubGeometry& ,
               const LocalCoefficients& ) {
    return LocalRefinedFunctionEvaluationHelper<isDGRefined>::
              evaluateValue(localFiniteElement, subElement, quadPos);
  }

  template<EvaluationType type, bool isDGRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::grad>* = nullptr>
  static std::vector<FieldVector<double,1> >
  evaluateRhs_(const LocalFiniteElement& localFiniteElement,
               unsigned int subElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& geometry,
               const SubGeometry& subGeometryInReferenceElement,
               const LocalCoefficients& localCoefficients) {
    return LocalRefinedFunctionEvaluationHelper<isDGRefined>::
              evaluateGrad(localFiniteElement,
                           subElement,
                           quadPos,
                           geometry,
                           subGeometryInReferenceElement,
                           localCoefficients.localSecondDirection()(quadPos));
  }
};

template<int dim, LinearIntegrationType integrationType>
struct LocalRefinedLinearTermFunctionEvaluation {
  constexpr static auto evaluationType = (integrationType ==
                      LinearIntegrationType::valueFunction)
                      ? EvaluationType::value : EvaluationType::grad;

  template<bool isDGRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients>
  static std::vector<FieldVector<double,1> >
  evaluate(const LocalFiniteElement& localFiniteElement,
           unsigned int subElement,
           const FieldVector<double, dim>& quadPos,
           const Geometry& geometry,
           const SubGeometry& subGeometryInReferenceElement,
           const LocalCoefficients& localCoefficients)
  {
    return evaluate_<evaluationType, isDGRefined>
            (localFiniteElement,
             subElement,
             quadPos,
             geometry,
             subGeometryInReferenceElement,
             localCoefficients);
  }

  private:

  template<EvaluationType type, bool isDGRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::value>* = nullptr>
  static std::vector<FieldVector<double,1> >
  evaluate_(const LocalFiniteElement& localFiniteElement,
            unsigned int subElement,
            const FieldVector<double, dim>& quadPos,
            const Geometry& ,
            const SubGeometry& ,
            const LocalCoefficients& ) {
    return LocalRefinedFunctionEvaluationHelper<isDGRefined>::
              evaluateValue(localFiniteElement, subElement, quadPos);
  }

  template<EvaluationType type, bool isDGRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::grad>* = nullptr>
  static std::vector<FieldVector<double,1> >
  evaluate_(const LocalFiniteElement& localFiniteElement,
            unsigned int subElement,
            const FieldVector<double, dim>& quadPos,
            const Geometry& geometry,
            const SubGeometry& subGeometryInReferenceElement,
            const LocalCoefficients& localCoefficients) {
    return LocalRefinedFunctionEvaluationHelper<isDGRefined>::
              evaluateGrad(localFiniteElement,
                           subElement,
                           quadPos,
                           geometry,
                           subGeometryInReferenceElement,
                           localCoefficients.localDirection()(quadPos));
  }
};


template<class BoundaryValue, class Index,
         typename std::enable_if<
                    std::is_arithmetic<std::decay_t<BoundaryValue>>::value>
                              ::type* = nullptr >
inline double evaluateBoundary(BoundaryValue boundary, Index)
{
  return boundary;
}

template<class BoundaryValue, class Index,
         typename std::enable_if<
                    is_vector<std::decay_t<BoundaryValue>>::value>
                              ::type* = nullptr >
inline double evaluateBoundary(const BoundaryValue& boundary, Index i)
{
  return boundary[i];
}

}}

#endif
