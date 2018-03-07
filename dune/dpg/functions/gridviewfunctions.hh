// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_DPG_GRIDVIEWFUNCTIONS_HH
#define DUNE_FUNCTIONS_DPG_GRIDVIEWFUNCTIONS_HH

#include <type_traits>

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

template<class Constant, class GridView>
ConstantGridViewFunction<Constant, GridView>
makeConstantGridViewFunction(Constant constant, const GridView&)
{
  return {constant};
}

}} // end namespace Dune::Functions

#endif
