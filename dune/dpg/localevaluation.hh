// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_LOCALEVALUATION_HH
#define DUNE_DPG_LOCALEVALUATION_HH

#include <algorithm>
#include <tuple>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include "assemble_types.hh"
#include "type_traits.hh"
#include <dune/common/std/type_traits.hh>

namespace Dune {
namespace detail {

/* We need to make this a class, as partial specializations of
 * function templates are not allowed. */
template<int dim, class LocalFiniteElement>
void evaluateLocalFunctionValue
                      (std::vector<FieldVector<double,1>>& values,
                       const LocalFiniteElement& localFiniteElement,
                       const FieldVector<double, dim>& quadPos)
{
  localFiniteElement.localBasis().evaluateFunction(quadPos, values);
}

template<int dim, class LocalFiniteElement, class Geometry>
void evaluateLocalFunctionGrad
                      (std::vector<FieldVector<double,1>>& derivatives,
                       const LocalFiniteElement& localFiniteElement,
                       const FieldVector<double, dim> & quadPos,
                       const Geometry& geometry,
                       const FieldVector<double, dim>& beta)
{
  const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);
  // The gradients of the shape functions on the reference element
  std::vector<FieldMatrix<double,1,dim>> referenceGradients;
  localFiniteElement.localBasis()
          .evaluateJacobian(quadPos, referenceGradients);

  // Compute the shape function gradients on the real element
  derivatives.resize(referenceGradients.size());
  std::transform(cbegin(referenceGradients),
                 cend(referenceGradients),
                 begin(derivatives),
                 [&](const FieldMatrix<double,1,dim>& referenceGradient)
                 {
                    FieldVector<double,dim> gradient;
                    jacobian.mv(referenceGradient[0], gradient);
                    return beta * gradient;
                 });
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
  static void
  evaluateLhs(std::vector<FieldVector<double,1>>& lhsValues,
              const LocalFiniteElement& localFiniteElement,
              const FieldVector<double, dim>& quadPos,
              const Geometry& geometry,
              const LocalCoefficients& localCoefficients)
  {
    evaluateLhs_<lhsType>
                (lhsValues,
                 localFiniteElement,
                 quadPos,
                 geometry,
                 localCoefficients);
  }

  template<class LocalFiniteElement, class Geometry, class LocalCoefficients>
  static void
  evaluateRhs(std::vector<FieldVector<double,1>>& rhsValues,
              const LocalFiniteElement& localFiniteElement,
              const FieldVector<double, dim>& quadPos,
              const Geometry& geometry,
              const LocalCoefficients& localCoefficients)
  {
    evaluateRhs_<rhsType>
                (rhsValues,
                 localFiniteElement,
                 quadPos,
                 geometry,
                 localCoefficients);
  }

  private:

  template<EvaluationType type,
           class LocalFiniteElement, class Geometry, class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::value>* = nullptr>
  static void
  evaluateLhs_(std::vector<FieldVector<double,1>>& lhsValues,
               const LocalFiniteElement& localFiniteElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& ,
               const LocalCoefficients& ) {
    evaluateLocalFunctionValue(lhsValues, localFiniteElement, quadPos);
  }

  template<EvaluationType type,
           class LocalFiniteElement, class Geometry, class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::grad>* = nullptr>
  static void
  evaluateLhs_(std::vector<FieldVector<double,1>>& lhsValues,
               const LocalFiniteElement& localFiniteElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& geometry,
               const LocalCoefficients& localCoefficients) {
    evaluateLocalFunctionGrad
             (lhsValues,
              localFiniteElement,
              quadPos,
              geometry,
              localCoefficients.localDirection()(quadPos));
  }

  template<EvaluationType type,
           class LocalFiniteElement, class Geometry, class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::value>* = nullptr>
  static void
  evaluateRhs_(std::vector<FieldVector<double,1>>& rhsValues,
               const LocalFiniteElement& localFiniteElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& ,
               const LocalCoefficients& ) {
    evaluateLocalFunctionValue(rhsValues, localFiniteElement, quadPos);
  }

  template<EvaluationType type,
           class LocalFiniteElement, class Geometry, class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::grad>* = nullptr>
  static void
  evaluateRhs_(std::vector<FieldVector<double,1>>& rhsValues,
               const LocalFiniteElement& localFiniteElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& geometry,
               const LocalCoefficients& localCoefficients) {
    evaluateLocalFunctionGrad
             (rhsValues,
              localFiniteElement,
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

  template<class LocalFiniteElement, class Geometry, class LocalCoefficients>
  static void
  evaluate(std::vector<FieldVector<double,1>>& values,
           const LocalFiniteElement& localFiniteElement,
           const FieldVector<double, dim>& quadPos,
           const Geometry& geometry,
           const LocalCoefficients& localCoefficients) {
    evaluate_<evaluationType>
             (values,
              localFiniteElement,
              quadPos,
              geometry,
              localCoefficients);
  }

  private:

  template<EvaluationType type,
           class LocalFiniteElement, class Geometry, class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::value>* = nullptr>
  static void
  evaluate_(std::vector<FieldVector<double,1>>& values,
            const LocalFiniteElement& localFiniteElement,
            const FieldVector<double, dim>& quadPos,
            const Geometry& ,
            const LocalCoefficients& ) {
    evaluateLocalFunctionValue(values, localFiniteElement, quadPos);
  }

  template<EvaluationType type,
           class LocalFiniteElement, class Geometry, class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::grad>* = nullptr>
  static void
  evaluate_(std::vector<FieldVector<double,1>>& values,
            const LocalFiniteElement& localFiniteElement,
            const FieldVector<double, dim>& quadPos,
            const Geometry& geometry,
            const LocalCoefficients& localCoefficients) {
    evaluateLocalFunctionGrad
             (values,
              localFiniteElement,
              quadPos,
              geometry,
              localCoefficients.localDirection()(quadPos));
  }
};

/* We need to make this a class, as partial specializations of
 * function templates are not allowed. */
template<bool isContinuouslyRefined>
struct LocalRefinedFunctionEvaluationHelper {

  template <int dim, class LocalFiniteElement,
            class Geometry, class SubGeometry>
  static void
  evaluateValue(std::vector<FieldVector<double,1>>& values,
                const LocalFiniteElement& localFiniteElement,
                unsigned int subElement,
                const FieldVector<double, dim>& quadPos);

  template <int dim, class LocalFiniteElement,
            class Geometry, class SubGeometry>
  static void
  evaluateGrad(std::vector<FieldVector<double,1>>& derivatives,
               const LocalFiniteElement& localFiniteElement,
               unsigned int subElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& geometry,
               const SubGeometry& subGeometryInReferenceElement,
               const FieldVector<double, dim>& beta);
};

template<>
struct LocalRefinedFunctionEvaluationHelper<false> {

  template<int dim, class LocalFiniteElement>
  static void
  evaluateValue(std::vector<FieldVector<double,1>>& values,
                const LocalFiniteElement& localFiniteElement,
                unsigned int,
                const FieldVector<double, dim>& quadPos)
  {
    localFiniteElement.localBasis().evaluateFunction(quadPos, values);
  }

  template<int dim, class LocalFiniteElement,
           class Geometry, class SubGeometry>
  static void
  evaluateGrad(std::vector<FieldVector<double,1>>& derivatives,
               const LocalFiniteElement& localFiniteElement,
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
    derivatives.resize(referenceGradients.size());
    std::transform(cbegin(referenceGradients),
                   cend(referenceGradients),
                   begin(derivatives),
                   [&](const FieldMatrix<double,1,dim>& referenceGradient)
                   {
                     FieldVector<double,dim> gradientRef, gradient;
                     jacobianSub.mv(referenceGradient[0], gradientRef);
                     jacobian.mv(gradientRef, gradient);
                     return beta * gradient;
                   });
  }
};

template<>
struct LocalRefinedFunctionEvaluationHelper<true> {

  template<int dim, class LocalFiniteElement>
  static void
  evaluateValue(std::vector<FieldVector<double,1>>& values,
                const LocalFiniteElement& localFiniteElement,
                unsigned int subElement,
                const FieldVector<double, dim>& quadPos)
  {
    localFiniteElement.localBasis().evaluateFunction(subElement, quadPos,
                                                     values);
  }

  template<int dim, class LocalFiniteElement,
           class Geometry, class SubGeometry>
  static void
  evaluateGrad(std::vector<FieldVector<double,1>>& derivatives,
               const LocalFiniteElement& localFiniteElement,
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
    derivatives.resize(referenceGradients.size());
    std::transform(cbegin(referenceGradients),
                   cend(referenceGradients),
                   begin(derivatives),
                   [&](const FieldMatrix<double,1,dim>& referenceGradient)
                   {
                     FieldVector<double,dim> gradientRef, gradient;
                     jacobianSub.mv(referenceGradient[0], gradientRef);
                     jacobian.mv(gradientRef, gradient);
                     return beta * gradient;
                   });
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

  template<bool isContinuouslyRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients>
  static void
  evaluateLhs(std::vector<FieldVector<double,1>>& lhsValues,
              const LocalFiniteElement& localFiniteElement,
              unsigned int subElement,
              const FieldVector<double, dim>& quadPos,
              const Geometry& geometry,
              const SubGeometry& subGeometryInReferenceElement,
              const LocalCoefficients& localCoefficients)
  {
    evaluateLhs_<lhsType, isContinuouslyRefined>
     (lhsValues,
      localFiniteElement,
      subElement,
      quadPos,
      geometry,
      subGeometryInReferenceElement,
      localCoefficients);
  }

  template<bool isContinuouslyRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients>
  static void
  evaluateRhs(std::vector<FieldVector<double,1>>& rhsValues,
              const LocalFiniteElement& localFiniteElement,
              unsigned int subElement,
              const FieldVector<double, dim>& quadPos,
              const Geometry& geometry,
              const SubGeometry& subGeometryInReferenceElement,
              const LocalCoefficients& localCoefficients)
  {
    evaluateRhs_<rhsType, isContinuouslyRefined>
     (rhsValues,
      localFiniteElement,
      subElement,
      quadPos,
      geometry,
      subGeometryInReferenceElement,
      localCoefficients);
  }

  private:

  template<EvaluationType type, bool isContinuouslyRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::value>* = nullptr>
  static void
  evaluateLhs_(std::vector<FieldVector<double,1>>& lhsValues,
               const LocalFiniteElement& localFiniteElement,
               unsigned int subElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& ,
               const SubGeometry& ,
               const LocalCoefficients& ) {
    LocalRefinedFunctionEvaluationHelper<isContinuouslyRefined>::
       evaluateValue(lhsValues, localFiniteElement, subElement, quadPos);
  }

  template<EvaluationType type, bool isContinuouslyRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::grad>* = nullptr>
  static void
  evaluateLhs_(std::vector<FieldVector<double,1>>& lhsValues,
               const LocalFiniteElement& localFiniteElement,
               unsigned int subElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& geometry,
               const SubGeometry& subGeometryInReferenceElement,
               const LocalCoefficients& localCoefficients) {
    LocalRefinedFunctionEvaluationHelper<isContinuouslyRefined>::
       evaluateGrad(lhsValues,
                    localFiniteElement,
                    subElement,
                    quadPos,
                    geometry,
                    subGeometryInReferenceElement,
                    localCoefficients.localDirection()(quadPos));
  }

  template<EvaluationType type, bool isContinuouslyRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::value>* = nullptr>
  static void
  evaluateRhs_(std::vector<FieldVector<double,1>>& rhsValues,
               const LocalFiniteElement& localFiniteElement,
               unsigned int subElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& ,
               const SubGeometry& ,
               const LocalCoefficients& ) {
    LocalRefinedFunctionEvaluationHelper<isContinuouslyRefined>::
       evaluateValue(rhsValues, localFiniteElement, subElement, quadPos);
  }

  template<EvaluationType type, bool isContinuouslyRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::grad>* = nullptr>
  static void
  evaluateRhs_(std::vector<FieldVector<double,1>>& rhsValues,
               const LocalFiniteElement& localFiniteElement,
               unsigned int subElement,
               const FieldVector<double, dim>& quadPos,
               const Geometry& geometry,
               const SubGeometry& subGeometryInReferenceElement,
               const LocalCoefficients& localCoefficients) {
    LocalRefinedFunctionEvaluationHelper<isContinuouslyRefined>::
       evaluateGrad(rhsValues,
                    localFiniteElement,
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

  template<bool isContinuouslyRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients>
  static void
  evaluate(std::vector<FieldVector<double,1>>& values,
           const LocalFiniteElement& localFiniteElement,
           unsigned int subElement,
           const FieldVector<double, dim>& quadPos,
           const Geometry& geometry,
           const SubGeometry& subGeometryInReferenceElement,
           const LocalCoefficients& localCoefficients)
  {
    evaluate_<evaluationType, isContinuouslyRefined>
     (values,
      localFiniteElement,
      subElement,
      quadPos,
      geometry,
      subGeometryInReferenceElement,
      localCoefficients);
  }

  private:

  template<EvaluationType type, bool isContinuouslyRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::value>* = nullptr>
  static void
  evaluate_(std::vector<FieldVector<double,1>>& values,
            const LocalFiniteElement& localFiniteElement,
            unsigned int subElement,
            const FieldVector<double, dim>& quadPos,
            const Geometry& ,
            const SubGeometry& ,
            const LocalCoefficients& ) {
    LocalRefinedFunctionEvaluationHelper<isContinuouslyRefined>::
       evaluateValue(values, localFiniteElement, subElement, quadPos);
  }

  template<EvaluationType type, bool isContinuouslyRefined,
           class LocalFiniteElement, class Geometry, class SubGeometry,
           class LocalCoefficients,
           std::enable_if_t<type == EvaluationType::grad>* = nullptr>
  static void
  evaluate_(std::vector<FieldVector<double,1>>& derivatives,
            const LocalFiniteElement& localFiniteElement,
            unsigned int subElement,
            const FieldVector<double, dim>& quadPos,
            const Geometry& geometry,
            const SubGeometry& subGeometryInReferenceElement,
            const LocalCoefficients& localCoefficients) {
    LocalRefinedFunctionEvaluationHelper<isContinuouslyRefined>::
       evaluateGrad(derivatives,
                    localFiniteElement,
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
