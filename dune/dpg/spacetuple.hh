// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_SPACETUPLE_HH
#define DUNE_DPG_SPACETUPLE_HH

#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/tupleutility.hh>

#include <dune/dpg/functions/concepts.hh>

namespace Dune {

  template<class... Spaces>
  using SpaceTuplePtr = std::shared_ptr<std::tuple<Spaces...>>;

  template<class Spaces, std::size_t offset, std::size_t count>
  class SpaceTupleView;

  template<class T>
  struct is_SpaceTuplePtr : std::false_type {};

  template<class... Spaces>
  struct is_SpaceTuplePtr<SpaceTuplePtr<Spaces...>> : std::true_type {};

  template<class... Spaces, std::size_t offset, std::size_t count>
  struct is_SpaceTuplePtr<SpaceTupleView<SpaceTuplePtr<Spaces...>,
                                         offset, count>>
      : is_SpaceTuplePtr<SpaceTuplePtr<Spaces...>> {};

  /**
   * Create a shared_ptr to a tuple of spaces over the same GridView
   */
  template<class... Spaces, class GridView>
  SpaceTuplePtr<Spaces...>
  make_space_tuple(const GridView& gridView)
  {
    static_assert(Concept::tupleEntriesModel<
        Functions::Concept::GeneralizedGlobalBasis<GridView>,
        std::tuple<Spaces...>>(),
        "Spaces need to model the GeneralizedGlobalBasis concept.");

    return std::make_shared<typename std::tuple<Spaces...>>(
        std::make_tuple(Spaces(gridView)...));
  }

  template<class Spaces, std::size_t offset, std::size_t count>
  class SpaceTupleView {
    using SpaceTuple = typename Spaces::element_type;
    static_assert(is_SpaceTuplePtr<Spaces>::value,
        "Spaces needs to be a SpaceTuplePtr!");
    static_assert(offset+count <= std::tuple_size<SpaceTuple>::value,
        "offset+count exceed size of underlying space tuple!");

    const Spaces spaces;

    SpaceTuple& getSpaceTuple() noexcept { return *spaces; }
    const SpaceTuple& getSpaceTuple() const noexcept { return *spaces; }

    template<std::size_t I, class Spaces_,
             std::size_t offset_, std::size_t count_>
    friend
    constexpr
    std::tuple_element_t<I, SpaceTupleView<Spaces_, offset_, count_>>&
    std::get(SpaceTupleView<Spaces_, offset_, count_>& spaces) noexcept;

    template<std::size_t I, class Spaces_,
             std::size_t offset_, std::size_t count_>
    friend
    constexpr
    std::tuple_element_t<I, SpaceTupleView<Spaces_, offset_, count_>>&&
    std::get(SpaceTupleView<Spaces_, offset_, count_>&& spaces) noexcept;

    template<std::size_t I, class Spaces_,
             std::size_t offset_, std::size_t count_>
    friend
    constexpr
    std::tuple_element_t<I, SpaceTupleView<Spaces_, offset_, count_>> const&
    std::get(SpaceTupleView<Spaces_, offset_, count_> const& spaces) noexcept;

    template<std::size_t I, class Spaces_,
             std::size_t offset_, std::size_t count_>
    friend
    constexpr
    std::tuple_element_t<I, SpaceTupleView<Spaces_, offset_, count_>> const&&
    std::get(SpaceTupleView<Spaces_, offset_, count_> const&& spaces) noexcept;

  public:

    using SpaceTuplePtr = Spaces;

    SpaceTupleView(const SpaceTuplePtr& spaces) : spaces(spaces) {}

    SpaceTuplePtr getSpaceTuplePtr() noexcept { return spaces; }

    // To make SpaceTupleView behave as a SpaceTuplePtr
    using element_type = SpaceTupleView;
    SpaceTupleView& operator*() noexcept
    { return *this; }
    const SpaceTupleView& operator*() const noexcept
    { return *this; }
    SpaceTupleView* operator->() noexcept
    { return this; }
    const SpaceTupleView* operator->() const noexcept
    { return this; }
  };

  // TODO: make sure that we do not end up with a matryoshka of SpaceTupleView's
  template<std::size_t offset, std::size_t count, class SpaceTuplePtr>
  SpaceTupleView<SpaceTuplePtr, offset, count>
  make_space_tuple_view(const SpaceTuplePtr& spaces)
  {
    return SpaceTupleView<SpaceTuplePtr, offset, count>(spaces);
  }

  namespace detail {
    template<template <class> class TE,
             class View, class Indices>
    struct SpaceTupleViewForEachImpl;

    template<template <class> class TE,
             class View, std::size_t... I>
    struct SpaceTupleViewForEachImpl<TE, View, std::index_sequence<I...>>
    {
      using Type
          = std::tuple<typename TE<std::tuple_element_t<I, View>>::Type...>;
    };
  }

} // end namespace Dune

namespace std {
  template<std::size_t I, class Spaces, std::size_t offset, std::size_t count>
  struct tuple_element<I, Dune::SpaceTupleView<Spaces, offset, count>>
      : tuple_element<I+offset, typename Spaces::element_type>
  {
    static_assert(
        I+offset < std::tuple_size<typename Spaces::element_type>::value,
        "Index I out of bounds!");
  };

  template<class Spaces, std::size_t offset, std::size_t count>
  struct tuple_size<Dune::SpaceTupleView<Spaces, offset, count>>
      : integral_constant<std::size_t, count> {};

  template<std::size_t I, class Spaces, std::size_t offset, std::size_t count>
  constexpr
  tuple_element_t<I, Dune::SpaceTupleView<Spaces, offset, count>>&
  get(Dune::SpaceTupleView<Spaces, offset, count>& spaces) noexcept
  {
    static_assert(
        I+offset < std::tuple_size<typename Spaces::element_type>::value,
        "Index I out of bounds!");
    return get<I+offset>(spaces.getSpaceTuple());
  }

  template<std::size_t I, class Spaces, std::size_t offset, std::size_t count>
  constexpr
  tuple_element_t<I, Dune::SpaceTupleView<Spaces, offset, count>>&&
  get(Dune::SpaceTupleView<Spaces, offset, count>&& spaces) noexcept
  {
    static_assert(
        I+offset < std::tuple_size<typename Spaces::element_type>::value,
        "Index I out of bounds!");
    return get<I+offset>(spaces.getSpaceTuple());
  }

  template<std::size_t I, class Spaces, std::size_t offset, std::size_t count>
  constexpr
  tuple_element_t<I, Dune::SpaceTupleView<Spaces, offset, count>> const&
  get(Dune::SpaceTupleView<Spaces, offset, count> const & spaces) noexcept
  {
    static_assert(
        I+offset < std::tuple_size<typename Spaces::element_type>::value,
        "Index I out of bounds!");
    return get<I+offset>(spaces.getSpaceTuple());
  }

  template<std::size_t I, class Spaces, std::size_t offset, std::size_t count>
  constexpr
  tuple_element_t<I, Dune::SpaceTupleView<Spaces, offset, count>> const&&
  get(Dune::SpaceTupleView<Spaces, offset, count> const && spaces) noexcept
  {
    static_assert(
        I+offset < std::tuple_size<typename Spaces::element_type>::value,
        "Index I out of bounds!");
    return get<I+offset>(spaces.getSpaceTuple());
  }
} // end namespace std

namespace Dune {
  template<template <class> class TE,
           std::size_t offset, std::size_t count, class SpaceTuplePtr>
  class ForEachType<TE, SpaceTupleView<SpaceTuplePtr, offset, count>>
  {
    using View = SpaceTupleView<SpaceTuplePtr, offset, count>;
    using Indices = std::make_index_sequence<count>;

  public:
    using Type
        = typename detail::SpaceTupleViewForEachImpl<TE, View, Indices>::Type;
  };

#ifndef DOXYGEN
  template<class Tuple, class Functor, std::size_t... I>
  inline auto genericTransformSpaceTupleBackendImpl(
      Tuple& t, Functor& f, const std::index_sequence<I...>& )
    -> std::tuple<decltype(f(std::get<I>(t)))...>
  {
    return std::tuple<decltype(f(std::get<I>(t)))...>(f(std::get<I>(t))...);
  }

  template<std::size_t offset, std::size_t count, class SpaceTuplePtr,
           class Functor>
  auto genericTransformSpaceTupleBackend(
      SpaceTupleView<SpaceTuplePtr, offset, count>& t, Functor& f)
  -> decltype(genericTransformSpaceTupleBackendImpl(t, f,
                                  std::make_index_sequence<count>{}))
  {
    return genericTransformSpaceTupleBackendImpl(t, f,
                                  std::make_index_sequence<count>{});
  }

  template<std::size_t offset, std::size_t count, class SpaceTuplePtr,
           class Functor>
  auto genericTransformSpaceTupleBackend(
      const SpaceTupleView<SpaceTuplePtr, offset, count>& t, Functor& f)
  -> decltype(genericTransformSpaceTupleBackendImpl(t, f,
                                  std::make_index_sequence<count>{}))
  {
    return genericTransformSpaceTupleBackendImpl(t, f,
                                  std::make_index_sequence<count>{});
  }

  template<class... Args, class Functor>
  auto genericTransformSpaceTupleBackend(std::tuple<Args...>& t, Functor& f) ->
    decltype(genericTransformSpaceTupleBackendImpl(t, f, std::index_sequence_for<Args...>{}))
  {
    return genericTransformSpaceTupleBackendImpl(t, f, std::index_sequence_for<Args...>{});
  }

  template<class... Args, class Functor>
  auto genericTransformSpaceTupleBackend(const std::tuple<Args...>& t, Functor& f) ->
    decltype(genericTransformSpaceTupleBackendImpl(t, f, std::index_sequence_for<Args...>{}))
  {
    return genericTransformSpaceTupleBackendImpl(t, f, std::index_sequence_for<Args...>{});
  }
#endif

  template<class SpaceTuple, class Functor>
  auto genericTransformSpaceTuple(SpaceTuple&& t, Functor&& f) ->
    decltype(genericTransformSpaceTupleBackend(t, f))
  {
    return genericTransformSpaceTupleBackend(t, f);
  }
} // End namespace Dune

#endif
