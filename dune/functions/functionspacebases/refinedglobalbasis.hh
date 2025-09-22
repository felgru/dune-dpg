// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDGLOBALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDGLOBALBASIS_HH

#include <type_traits>

#include <dune/common/concept.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/version.hh>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/concepts.hh>
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,11)
#include <dune/functions/functionspacebases/containerdescriptors.hh>
#endif
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/refinedlocalview.hh>



namespace Dune {
namespace Functions {



/**
 * \brief Global basis for given prebasis
 *
 * This class implements the interface of a global basis
 * using the details from a given prebasis. Hence
 * it serves as an example for this interface.
 *
 * If you want to implement your own global basis, it may be
 * better to implement a prebasis instead. On the one hand
 * this needs less boiler-plate code. On the other hand
 * it makes your implementation composable and thus much
 * more flexible. That is, you can reuse your prebasis
 * as one part in a larger product space by plugging it
 * e.g. into a CompositePreBasis of PowerPreBasis.
 * The actual global basis for your FooPreBasis is
 * then obtained by using RefinedGlobalBasis<FooPreBasis>.
 *
 * \tparam PB  Prebasis providing the implementation details
 */
template<class PB>
class RefinedGlobalBasis
{
public:

  //! Prebasis providing the implementation details
  using PreBasis = PB;

  //! The grid view that the FE space is defined on
  using GridView = typename PreBasis::GridView;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = std::conditional_t<
      (PreBasis::multiIndexBufferSize == 1),
      FlatMultiIndex<std::size_t>,
      Dune::ReservedVector<std::size_t, PreBasis::multiIndexBufferSize>>;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Type of the local view on the restriction of the basis to a single element
  using LocalView = RefinedLocalView<RefinedGlobalBasis<PreBasis>>;

  //! Type used for prefixes handed to the size() method
  using SizePrefix
    = Dune::ReservedVector<std::size_t, PreBasis::multiIndexBufferSize>;

  /**
   * \brief Constructor
   *
   * \tparam T Argument list for PreBasis
   * \param t Argument list for PreBasis
   *
   * This will forward all arguments to the constructor of PreBasis
   */
  template<class... T,
    disableCopyMove<RefinedGlobalBasis, T...> = 0,
    enableIfConstructible<PreBasis, T...> = 0>
  RefinedGlobalBasis(T&&... t) :
    preBasis_(std::forward<T>(t)...)
  {
    static_assert(models<Concept::PreBasis<GridView>, PreBasis>(), "Type passed to RefinedGlobalBasis does not model the PreBasis concept.");
    preBasis_.initializeIndices();
  }

  /**
   * \brief Constructor from a PreBasis factory
   *
   * \param gridView  The GridView this basis is based on
   * \param factory  A factory functor that gets the `gridView` and returns a `PreBasis`
   */
  template<class PreBasisFactory>
  RefinedGlobalBasis(const GridView& gridView, PreBasisFactory&& factory) :
    preBasis_(factory(gridView))
  {
    static_assert(models<Concept::PreBasis<GridView>, PreBasis>(), "Type passed to RefinedGlobalBasis does not model the PreBasis concept.");
    preBasis_.initializeIndices();
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return preBasis_.gridView();
  }

  //! Obtain the prebasis providing the implementation details
  const PreBasis& preBasis() const
  {
    return preBasis_;
  }

  /**
   * \brief Update the stored grid view
   *
   * This will update the indexing information of the global basis.
   * It must be called if the grid has changed.
   */
  void update(const GridView & gv)
  {
    preBasis_.update(gv);
    preBasis_.initializeIndices();
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return preBasis_.dimension();
  }

  //! Return number of possible values for next position in empty multi index
  size_type size() const
  {
    return preBasis_.size();
  }

  //! Return number of possible values for next position in multi index
  size_type size(const SizePrefix& prefix) const
  {
    return preBasis_.size(prefix);
  }

  //! Return local view for basis
  LocalView localView() const
  {
    return LocalView(*this);
  }

  //! Return *this because we are not embedded in a larger basis
  const RefinedGlobalBasis& rootBasis() const
  {
    return *this;
  }

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,11)
  //! Return the associated container descriptor
  auto containerDescriptor() const
  {
    if constexpr (requires(PreBasis pb){ pb.containerDescriptor(); })
      return preBasis_.containerDescriptor();
    else
      return ContainerDescriptors::Unknown{};
  }
#endif

protected:
  PreBasis preBasis_;
};



template<class PreBasis>
RefinedGlobalBasis(PreBasis&&) -> RefinedGlobalBasis<std::decay_t<PreBasis>>;

template<class GridView, class PreBasisFactory>
RefinedGlobalBasis(const GridView& gv, PreBasisFactory&& f) -> RefinedGlobalBasis<std::decay_t<decltype(f(gv))>>;



namespace BasisFactory {

template<class GridView, class PreBasisFactory>
auto makeRefiendBasis(const GridView& gridView,
                      PreBasisFactory&& preBasisFactory)
{
  return RefinedGlobalBasis(preBasisFactory(gridView));
}

} // end namespace BasisFactory



} // end namespace Functions
} // end namespace Dune



#endif
