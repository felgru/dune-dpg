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
template<int dim, EvaluationType type>
struct LocalFunctionEvaluation {

  template <class LocalFiniteElement, class Geometry>
  std::vector<FieldVector<double,1> >
  operator() (const LocalFiniteElement& localFiniteElement,
              const FieldVector<double, dim>& quadPos,
              const Geometry& geometry,
              const FieldVector<double, dim>& beta) const;
};

template<int dim>
struct LocalFunctionEvaluation<dim, EvaluationType::value> {

  template <class LocalFiniteElement, class Geometry>
  std::vector<FieldVector<double,1> > operator()
                      (const LocalFiniteElement& localFiniteElement,
                       const FieldVector<double, dim>& quadPos,
                       const Geometry&,
                       const FieldVector<double, dim>&) const
  {
    // values of the shape functions
    std::vector<FieldVector<double,1> > values;
    localFiniteElement.localBasis().evaluateFunction(quadPos, values);
    return values;
  }
};

template<int dim>
struct LocalFunctionEvaluation<dim, EvaluationType::grad> {

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
         bool isDGRefined>
struct LocalRefinedFunctionEvaluation {

  template <class LocalFiniteElement, class Geometry, class SubGeometry>
  std::vector<FieldVector<double,1> >
  operator() (const LocalFiniteElement& localFiniteElement,
              unsigned int subElement,
              const FieldVector<double, dim>& quadPos,
              const Geometry& geometry,
              const SubGeometry& subGeometryInReferenceElement,
              const FieldVector<double, dim>& beta) const;
};

template<int dim>
struct LocalRefinedFunctionEvaluation<dim, EvaluationType::value, false> {

  template <class LocalFiniteElement, class Geometry, class SubGeometry>
  std::vector<FieldVector<double,1> > operator()
                      (const LocalFiniteElement& localFiniteElement,
                       unsigned int,
                       const FieldVector<double, dim>& quadPos,
                       const Geometry&,
                       const SubGeometry&,
                       const FieldVector<double, dim>&) const
  {
    // values of the shape functions
    std::vector<FieldVector<double,1> > values;
    localFiniteElement.localBasis().evaluateFunction(quadPos, values);
    return values;
  }
};

template<int dim>
struct LocalRefinedFunctionEvaluation<dim, EvaluationType::grad, false> {

  template <class LocalFiniteElement, class Geometry, class SubGeometry>
  std::vector<FieldVector<double,1> > operator()
                      (const LocalFiniteElement& localFiniteElement,
                       unsigned int,
                       const FieldVector<double, dim> & quadPos,
                       const Geometry& geometry,
                       const SubGeometry& subGeometryInReferenceElement,
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

template<int dim>
struct LocalRefinedFunctionEvaluation<dim, EvaluationType::value, true> {

  template <class LocalFiniteElement, class Geometry, class SubGeometry>
  std::vector<FieldVector<double,1> > operator()
                      (const LocalFiniteElement& localFiniteElement,
                       unsigned int subElement,
                       const FieldVector<double, dim>& quadPos,
                       const Geometry& geometry,
                       const SubGeometry& subGeometryInReferenceElement,
                       const FieldVector<double, dim>&) const
  {
    // values of the shape functions
    std::vector<FieldVector<double,1> > values;
    localFiniteElement.localBasis().evaluateFunction(subElement, quadPos, values);
    return values;
  }
};

template<int dim>
struct LocalRefinedFunctionEvaluation<dim, EvaluationType::grad, true> {

  template <class LocalFiniteElement, class Geometry, class SubGeometry>
  std::vector<FieldVector<double,1> > operator()
                      (const LocalFiniteElement& localFiniteElement,
                       unsigned int subElement,
                       const FieldVector<double, dim> & quadPos,
                       const Geometry& geometry,
                       const SubGeometry& subGeometryInReferenceElement,
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
