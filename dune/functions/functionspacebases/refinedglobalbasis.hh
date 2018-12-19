// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDGLOBALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDGLOBALBASIS_HH

#include <dune/common/concept.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/version.hh>

#include <dune/functions/common/type_traits.hh>
#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
#include <dune/functions/functionspacebases/defaultlocalindexset.hh>
#endif
#include <dune/functions/functionspacebases/refinedlocalview.hh>
#include <dune/functions/functionspacebases/concepts.hh>



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

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
  //! The empty prefix path that identifies the root in the local ansatz tree
  using PrefixPath = TypeTree::HybridTreePath<>;
#endif

  //! The grid view that the FE space is defined on
  using GridView = typename PreBasis::GridView;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = typename PreBasis::MultiIndex;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Type of the local view on the restriction of the basis to a single element
  using LocalView = RefinedLocalView<RefinedGlobalBasis<PreBasis>>;

  //! Node index set provided by PreBasis
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,7)
  using NodeIndexSet = typename PreBasis::IndexSet;
#else
  using NodeIndexSet = typename PreBasis::template IndexSet<PrefixPath>;
#endif

  //! Type used for prefixes handed to the size() method
  using SizePrefix = typename PreBasis::SizePrefix;

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
  //! Type of local index set exported by localIndexSet()
  using LocalIndexSet = DefaultLocalIndexSet<LocalView, NodeIndexSet>;
#endif


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
#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
    , prefixPath_()
#endif
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

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
  //! Return local index set for basis
  LocalIndexSet localIndexSet() const
  {
    return LocalIndexSet(preBasis_.template indexSet<PrefixPath>());
  }
#endif

  //! Return *this because we are not embedded in a larger basis
  const RefinedGlobalBasis& rootBasis() const
  {
    return *this;
  }

#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
  //! Return empty path, because this is the root in the local ansatz tree
  const PrefixPath& prefixPath() const
  {
    return prefixPath_;
  }
#endif

protected:
  PreBasis preBasis_;
#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,7)
  PrefixPath prefixPath_;
#endif
};



} // end namespace Functions
} // end namespace Dune



#endif
