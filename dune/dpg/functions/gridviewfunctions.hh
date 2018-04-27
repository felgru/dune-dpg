// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_DPG_GRIDVIEWFUNCTIONS_HH
#define DUNE_FUNCTIONS_DPG_GRIDVIEWFUNCTIONS_HH

#include <type_traits>

#include <dune/functions/common/functionconcepts.hh>
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
   * \brief Obtain local contex this LocalFunction is bound to
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
   * \brief Obtain local contex this LocalFunction is bound to
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

template<class Factor, class GridView,
         typename std::enable_if<
                  std::is_arithmetic<Factor>::value>
                            ::type* = nullptr >
ConstantGridViewFunction<Factor, GridView>
make_GridViewFunction(Factor factor, const GridView&)
{
  return {factor};
}

template<class Factor, class GridView,
         class std::enable_if<Concept::isGridViewFunction<Factor,
                      double(typename GridView::template Codim<0>
                                     ::Geometry::GlobalCoordinate),
                      GridView>()>
                  ::type* = nullptr>
Factor make_GridViewFunction(Factor factor, const GridView&)
{
  return factor;
}

}} // end namespace Dune::Functions

#endif
