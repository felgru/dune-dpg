#ifndef DUNE_DPG_FUNCTIONS_NORMEDSPACES_HH
#define DUNE_DPG_FUNCTIONS_NORMEDSPACES_HH

#include <memory>
#include <tuple>
#include <type_traits>

#include <dune/dpg/type_traits.hh>

#include <dune/functions/functionspacebases/normalizedbasisadaptor.hh>
#include <dune/functions/functionspacebases/normalizedrefinedbasisadaptor.hh>

namespace Dune {

/**
 * Create a shared_ptr to a tuple of normed spaces
 */
template<class InnerProduct>
auto make_normalized_space_tuple(const InnerProduct& innerProduct)
{
  using WrappedSpaces = typename InnerProduct::TestSpaces;

  static_assert(std::tuple_size<WrappedSpaces>::value == 1,
      "make_normalized_space_tuple only implemented for tuples with one space.");
  using NormedSpace = std::conditional_t<
    is_RefinedFiniteElement<std::tuple_element_t<0, WrappedSpaces>>::value,
    Functions::NormalizedRefinedBasis<InnerProduct>,
    Functions::NormalizedBasis<InnerProduct>>;

  return std::make_shared<std::tuple<NormedSpace>>(
      std::make_tuple(NormedSpace(innerProduct)));
}

} // end namespace Dune

#endif
