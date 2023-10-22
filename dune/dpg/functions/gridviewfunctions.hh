// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_DPG_GRIDVIEWFUNCTIONS_HH
#define DUNE_FUNCTIONS_DPG_GRIDVIEWFUNCTIONS_HH

#include <type_traits>
#include <utility>

#include <dune/common/fvector.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/common/typeutilities.hh>
#include <dune/dpg/quadratureorder.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>

namespace Dune {
namespace Functions {

template<class Range, class LocalDomain, class LocalContext>
class LocalConstantGridViewFunction
{
  using Domain = LocalDomain;

public:
  LocalConstantGridViewFunction(Range constant) : constant(constant) {}

  Range operator() (const Domain&) const
  {
    return constant;
  }

  void bind(const LocalContext& localContext)
  {
    context = &localContext;
  }

  /**
   * \brief Unbind from local context
   */
  void unbind()
  {
    context = nullptr;
  }

  /**
   * \brief Obtain the local context this LocalFunction is bound to
   */
  const LocalContext& localContext() const
  {
    return *context;
  }

private:
  const Range constant;
  const LocalContext* context = nullptr;
};

template<class Range, class GridView>
class ConstantGridViewFunction
{
  using EntitySet = GridViewEntitySet<GridView, 0>;

  using LocalCoordinate = typename EntitySet::LocalCoordinate;

  // using LocalSignature = typename Range(LocalCoordinate);

public:
  using Element = typename EntitySet::Element;

  using LocalFunction
    = LocalConstantGridViewFunction<Range, LocalCoordinate, Element>;

  ConstantGridViewFunction(Range constant) : constant(constant) {}

  friend LocalFunction
    localFunction(const ConstantGridViewFunction& t)
  {
    return LocalFunction(t.constant);
  }

private:
  const Range constant;
};

template<class Constant, class GridView>
ConstantGridViewFunction<std::decay_t<Constant>, GridView>
makeConstantGridViewFunction(Constant&& constant, const GridView&)
{
  return {std::forward<Constant>(constant)};
}

namespace detail {
  template<class T>
  using has_localFunction_t
      = decltype(localFunction(std::declval<const T&>()));

  template<class GridView, class T,
           typename std::enable_if<Std::is_detected_v<has_localFunction_t, T>>
                        ::type* = nullptr>
  T toGridViewFunction(T&& t) {
    return std::forward<T>(t);
  }

  template<class GridView, class T,
           typename std::enable_if<
                    std::is_arithmetic<std::decay_t<T>>::value>
                              ::type* = nullptr >
  ConstantGridViewFunction<std::decay_t<T>, GridView>
  toGridViewFunction(T&& t) {
    return {std::forward<T>(t)};
  }

  template<class GridView, class T, int dim,
           typename std::enable_if<
                    std::is_arithmetic<T>::value>
                              ::type* = nullptr >
  ConstantGridViewFunction<FieldVector<T,dim>, GridView>
  toGridViewFunction(FieldVector<T,dim>&& t) {
    return {std::move(t)};
  }

  template<class GridView, class T, int dim,
           typename std::enable_if<
                    std::is_arithmetic<T>::value>
                              ::type* = nullptr >
  ConstantGridViewFunction<FieldVector<T,dim>, GridView>
  toGridViewFunction(const FieldVector<T,dim>& t) {
    return {t};
  }
}

template<class Range, class LocalDomain, class LocalContext>
class LocalZeroGridViewFunction
{
  using Domain = LocalDomain;

public:
  LocalZeroGridViewFunction() {}

  Range operator() (const Domain&) const
  {
    return Range(0.);
  }

  //! Create a derivative grid-function returning constant 0.
  friend LocalZeroGridViewFunction<Range, LocalDomain, LocalContext>
    derivative(const LocalZeroGridViewFunction&)
  {
    return {};
  }

  void bind(const LocalContext& localContext)
  {
    context = &localContext;
  }

  /**
   * \brief Unbind from local context
   */
  void unbind()
  {
    context = nullptr;
  }

  /**
   * \brief Obtain the local context this LocalFunction is bound to
   */
  const LocalContext& localContext() const
  {
    return *context;
  }

private:
  const LocalContext* context = nullptr;
};

template<class Range, class GridView>
class ZeroGridViewFunction
{
  using EntitySet = GridViewEntitySet<GridView, 0>;
  using Domain
      = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
  using LocalCoordinate = typename EntitySet::LocalCoordinate;
  // using LocalSignature = typename Range(LocalCoordinate);

public:
  using Element = typename EntitySet::Element;

  using LocalFunction
    = LocalZeroGridViewFunction<Range, LocalCoordinate, Element>;

  ZeroGridViewFunction(const GridView& gridView) :
    entitySet_(gridView)
  {}

  Range operator()(const Domain& x) const
  {
    return Range(0.);
  }

  friend LocalFunction localFunction(const ZeroGridViewFunction&)
  {
    return LocalFunction();
  }

  //! Create a derivative grid-function returning constant 0.
  friend ZeroGridViewFunction<Range, GridView>
    derivative(const ZeroGridViewFunction& t)
  {
    return {t.entitySet_.gridView()};
  }

  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

private:
  EntitySet entitySet_;
};

template<class Range, class LocalDomain, class LocalContext, class F>
class LocalPiecewiseConstantGridViewFunction
{
  using Domain = LocalDomain;

public:
  LocalPiecewiseConstantGridViewFunction(Range constant) : constant(constant) {}

  template<class FT,
           disableCopyMove<LocalPiecewiseConstantGridViewFunction, FT> = 0>
  LocalPiecewiseConstantGridViewFunction(FT&& f) :
    f_(std::forward<FT>(f))
  {}

  Range operator() (const Domain&) const
  {
    return constant;
  }

  //! Create a derivative grid-function returning constant 0.
  friend LocalZeroGridViewFunction<Range, LocalDomain, LocalContext>
    derivative(const LocalPiecewiseConstantGridViewFunction&)
  {
    return {};
  }

  void bind(const LocalContext& localContext)
  {
    context = &localContext;
    const Domain centerOfLocalContext = localContext.geometry().center();
    constant = f_(centerOfLocalContext);
  }

  /**
   * \brief Unbind from local context
   */
  void unbind()
  {
    context = nullptr;
  }

  /**
   * \brief Obtain the local context this LocalFunction is bound to
   */
  const LocalContext& localContext() const
  {
    return *context;
  }

private:
  F f_;
  Range constant;
  const LocalContext* context = nullptr;
};

template<class Range, class GridView, class F>
class PiecewiseConstantGridViewFunction
{
  using EntitySet = GridViewEntitySet<GridView, 0>;
  using Domain
      = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
  using LocalCoordinate = typename EntitySet::LocalCoordinate;
  // using LocalSignature = typename Range(LocalCoordinate);

public:
  using Element = typename EntitySet::Element;

  using LocalFunction
    = LocalPiecewiseConstantGridViewFunction<Range, LocalCoordinate,
                                             Element, F>;

  template<class FT>
  PiecewiseConstantGridViewFunction(FT&& f, const GridView& gridView) :
    f_(std::forward<FT>(f)),
    entitySet_(gridView)
  {}

  Range operator()(const Domain& x) const
  {
    return f_(x);
  }

  friend LocalFunction
    localFunction(const PiecewiseConstantGridViewFunction& t)
  {
    return LocalFunction(t.f_);
  }

  //! Create a derivative grid-function returning constant 0.
  friend ZeroGridViewFunction<Range, GridView>
    derivative(const PiecewiseConstantGridViewFunction& t)
  {
    return {t.entitySet_.gridView()};
  }

  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

private:
  F f_;
  EntitySet entitySet_;
};

template<class F, class GridView>
PiecewiseConstantGridViewFunction<
  typename std::result_of<F(typename GridView::template Codim<0>::Geometry::GlobalCoordinate)>::type,  // Range
  GridView,
  typename std::decay<F>::type >                                                                      // Raw type of F (without & or &&)
  makePiecewiseConstantGridViewFunction(F&& f, const GridView& gridView)
{
  using Domain = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
  using Range = typename std::result_of<F(Domain)>::type;
  using FRaw = typename std::decay<F>::type;

  return PiecewiseConstantGridViewFunction<Range, GridView, FRaw>
                                          (std::forward<F>(f), gridView);
}

} // end namespace Dune::Functions

template<class Range, class LocalDomain, class LocalContext>
struct requiredQuadratureOrder<Functions::LocalConstantGridViewFunction
                               <Range, LocalDomain, LocalContext>>
  : std::integral_constant<unsigned int, 0> {};

template<class Range, class LocalDomain, class LocalContext, class F>
struct requiredQuadratureOrder<Functions::LocalPiecewiseConstantGridViewFunction
                               <Range, LocalDomain, LocalContext, F>>
  : std::integral_constant<unsigned int, 0> {};

} // end namespace Dune

#endif
