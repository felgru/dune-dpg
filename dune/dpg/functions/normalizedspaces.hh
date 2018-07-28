#ifndef DUNE_DPG_FUNCTIONS_NORMALIZEDSPACES_HH
#define DUNE_DPG_FUNCTIONS_NORMALIZEDSPACES_HH

#include <memory>
#include <tuple>
#include <type_traits>

#include <dune/dpg/type_traits.hh>

#include <dune/functions/functionspacebases/normalizedbasisadaptor.hh>
#include <dune/functions/functionspacebases/normalizedrefinedbasisadaptor.hh>

namespace Dune {

/**
 * Create a normalized space
 */
template<class InnerProduct>
auto make_normalized_space(const InnerProduct& innerProduct)
{
  using WrappedSpaces = typename InnerProduct::TestSpaces;

  static_assert(std::tuple_size<WrappedSpaces>::value == 1,
      "make_normalized_space_tuple only implemented for tuples with one space.");
  using NormalizedSpace = std::conditional_t<
    is_RefinedFiniteElement<std::tuple_element_t<0, WrappedSpaces>>::value,
    Functions::NormalizedRefinedBasis<InnerProduct>,
    Functions::NormalizedBasis<InnerProduct>>;

  return NormalizedSpace(innerProduct);
}

/**
 * Create a shared_ptr to a tuple of normalized spaces
 */
template<class InnerProduct>
auto make_normalized_space_tuple(const InnerProduct& innerProduct)
{
  auto normalizedSpace = make_normalized_space(innerProduct);
  using NormalizedSpace = decltype(normalizedSpace);

  return std::make_shared<std::tuple<NormalizedSpace>>(
      std::make_tuple(std::move(normalizedSpace)));
}

} // end namespace Dune

#endif
