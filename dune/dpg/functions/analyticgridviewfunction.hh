// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_DPG_GRIDVIEWFUNCTION_HH
#define DUNE_FUNCTIONS_DPG_GRIDVIEWFUNCTION_HH

#include <dune/dpg/quadratureorder.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

namespace Dune {
namespace Functions {

namespace Imp {

template<unsigned int quadratureOrder, class Signature, class GV, class FLocal, template<class> class DerivativeTraits=DefaultDerivativeTraits>
class LocalAnalyticGridViewFunctionWithQuadratureOrder;

template<unsigned int quadratureOrder, class Range, class LocalDomain, class GV, class F, template<class> class DerivativeTraits>
class LocalAnalyticGridViewFunctionWithQuadratureOrder<quadratureOrder, Range(LocalDomain), GV, F, DerivativeTraits>
  : public LocalAnalyticGridViewFunction<Range(LocalDomain), GV, F, DerivativeTraits>
{
public:

  template<class FT, disableCopyMove<LocalAnalyticGridViewFunctionWithQuadratureOrder, FT> = 0>
  LocalAnalyticGridViewFunctionWithQuadratureOrder(FT&& f) :
    LocalAnalyticGridViewFunction
      <Range(LocalDomain), GV, F, DerivativeTraits>(std::forward<FT>(f))
  {}
};

} // end namespace Imp

template<unsigned int quadratureOrder, class Signature, class GV, class F, template<class> class DerivativeTraits=DefaultDerivativeTraits>
class AnalyticGridViewFunctionWithQuadratureOrder;

/**
 * \brief Class wrapping any differentiable function as grid function
 *
 * \ingroup FunctionImplementations
 */
template<unsigned int quadratureOrder, class Range, class Domain, class GV, class F, template<class> class DerivativeTraits>
class AnalyticGridViewFunctionWithQuadratureOrder<quadratureOrder, Range(Domain), GV, F, DerivativeTraits>
{
public:
  using Signature = Range(Domain);
  using RawSignature = typename SignatureTraits<Signature>::RawSignature;
  using DerivativeSignature = typename DerivativeTraits<RawSignature>::Range(Domain);

  using GridView = GV;
  using EntitySet = GridViewEntitySet<GridView, 0>;
  using Element = typename EntitySet::Element;
  using Geometry = typename Element::Geometry;

  // Use the indirection via derivativeIfImplemented to also support
  // function types F that do not implement derivative. In this case
  // the interface type DifferentiableFunction is using a dummy for
  // the derivative type.
  using DerivativeDummy = DifferentiableFunction<DerivativeSignature>;
  using GlobalRawDerivative = decltype(Imp::derivativeIfImplemented<DerivativeDummy, F>(std::declval<F>()));
  using Derivative = AnalyticGridViewFunction<DerivativeSignature, GridView, GlobalRawDerivative, DerivativeTraits>;

  using LocalDomain = typename EntitySet::LocalCoordinate;
  using LocalFunction = typename Imp::LocalAnalyticGridViewFunctionWithQuadratureOrder<quadratureOrder, Range(LocalDomain), GridView, F, LocalDerivativeTraits<EntitySet, DerivativeTraits>::template Traits>;

  template<class FT>
  AnalyticGridViewFunctionWithQuadratureOrder
  (FT&& f, const GridView& gridView) :
    f_(std::forward<FT>(f)),
    entitySet_(gridView)
  {}

  Range operator()(const Domain& x) const
  {
    return f_(x);
  }

  friend Derivative derivative
    (const AnalyticGridViewFunctionWithQuadratureOrder& t)
  {
    return Derivative(Imp::derivativeIfImplemented<DerivativeDummy, F>(t.f_), t.entitySet_.gridView());
  }

  friend LocalFunction localFunction
    (const AnalyticGridViewFunctionWithQuadratureOrder& t)
  {
    return LocalFunction(t.f_);
  }

  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

private:
  F f_;
  EntitySet entitySet_;
};

template<unsigned int quadratureOrder, class F, class GridView>
AnalyticGridViewFunctionWithQuadratureOrder<quadratureOrder,
  typename std::result_of<F(typename GridView::template Codim<0>::Geometry::GlobalCoordinate)>::type  // Range
  (typename GridView::template Codim<0>::Geometry::GlobalCoordinate),                                 // Domain
  GridView,
  typename std::decay<F>::type >                                                                      // Raw type of F (without & or &&)
  makeAnalyticGridViewFunctionWithQuadratureOrder(F&& f,
                                                  const GridView& gridView)
{
  using Domain = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
  using Range = typename std::result_of<F(Domain)>::type;
  using FRaw = typename std::decay<F>::type;

  return AnalyticGridViewFunctionWithQuadratureOrder<quadratureOrder,
            Range(Domain), GridView, FRaw>(std::forward<F>(f), gridView);
}

} // end namespace Functions

template<unsigned int quadratureOrder, class Range, class LocalDomain, class GV, class F, template<class> class DerivativeTraits>
struct requiredQuadratureOrder<Functions::Imp::LocalAnalyticGridViewFunctionWithQuadratureOrder<quadratureOrder, Range(LocalDomain), GV, F, DerivativeTraits>>
  : std::integral_constant<unsigned int, quadratureOrder> {};

} // end namespace Dune

#endif
