// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_LINEARFORM_HH
#define DUNE_DPG_LINEARFORM_HH

#include <tuple>
#include <vector>
#include <memory>
#include <type_traits>
#include <cassert>

#include <dune/common/hybridutilities.hh>

#include "assemble_helper.hh"
#include "assemble_types.hh"
#include "linearintegralterm.hh"
#include "spacetuple.hh"
#include "type_traits.hh"

namespace Dune {

  /**
   * \brief This class describes a linear form.
   *
   * \tparam Sps             shared_ptr of a tuple of (test) spaces
   * \tparam LinearTerms     tuple of LinearIntegralTerm
   */
  template<class Sps, class LinearTerms>
  class LinearForm
  {
    static_assert(is_SpaceTuplePtr<Sps>::value,
        "Sps needs to be a SpaceTuplePtr!");
  public:
    //! shared pointer to tuple type of test spaces
    using SpacesPtr = Sps;
    //! tuple type of test spaces
    using Spaces = typename SpacesPtr::element_type;
    static_assert(Concept::tupleEntriesModel<
        Functions::Concept::GeneralizedGlobalBasis<
            typename std::tuple_element_t<0, Spaces>::GridView>,
        Spaces>(),
        "Spaces need to model the GeneralizedGlobalBasis concept.");
    //! tuple type for the local views of the test spaces
    using LocalViews = detail::getLocalViews_t<Spaces>;

    LinearForm () = delete;
    /**
     * \brief constructor for LinearForm
     *
     * \note For your convenience, use make_LinearForm() instead.
     */
    constexpr LinearForm (const SpacesPtr&    spaces,
                          const LinearTerms&  terms)
               : spaces(spaces),
                 terms(terms),
                 localViews(nullptr)
    { }

    /**
     * \brief Compute the rhs vector for a single element.
     *
     * \pre bind() has to be called on the current element before.
     *
     * \param[out] elementVector will hold the local vector
     */
    template <class VectorType>
    void getLocalVector(VectorType& elementVector) const
    {
      // Set all entries to zero
      elementVector.resize(localTotalSpaceSize);
      elementVector = 0;

      Hybrid::forEach(terms,
              detail::getLocalVectorHelper
                      <VectorType,
                       Spaces>
                      (*localViews,
                       elementVector,
                       localSpaceOffsets));
    }

    /**
     * \brief Bind the LinearForm to the given localViews.
     *
     * \pre The given localViews have to be bound to the same element.
     */
    void bind(LocalViews& lv)
    {
      constexpr bool usesOptimalTestBasis =
            is_OptimalTestSpace<
                typename std::tuple_element<std::tuple_size<Spaces>::value-1,
                                            Spaces>::type
            >::value;

      localViews = std::addressof(lv);

      const auto& e = std::get<0>(lv).element();
      Hybrid::forEach(terms, [&](auto& t) { t.term.bind(e); });

      constexpr size_t size = std::tuple_size<LocalViews>::value;

      /* set up local offsets */
      if(!usesOptimalTestBasis) {
        localTotalSpaceSize = detail::computeOffsets(localSpaceOffsets, lv);
      } else { /* DPG formulation */
        for(size_t i=0; i<size; ++i)
        {
          localSpaceOffsets[i] = 0;
        }
        localTotalSpaceSize = std::get<size-1>(lv).size();
      }
    }

    /**
     * \brief Does exactly what it says on the tin.
     */
    SpacesPtr getSpaces() const
    { return spaces; }

    /**
     * \brief Does exactly what it says on the tin.
     */
    const LinearTerms& getTerms() const
    {
      return terms;
    }

    using SpaceIndexArray = size_t[std::tuple_size<Spaces>::value];

    /**
     * \brief Does exactly what it says on the tin.
     */
    const SpaceIndexArray& getLocalSpaceOffsets() const
    { return localSpaceOffsets; }

  private:
    SpacesPtr    spaces;
    LinearTerms  terms;

    size_t localSpaceOffsets[std::tuple_size<Spaces>::value];
    size_t localTotalSpaceSize;

    LocalViews* localViews;
  };

/**
 * \brief Creates a LinearForm (i.e. for a saddlepoint formulation),
 *        deducing the target type from the types of arguments.
 *
 * \param spaces          a shared_ptr to a tuple of spaces
 * \param terms           a tuple of LinearIntegralTerm
 */
template<class SpacesPtr, class LinearTerms>
auto make_LinearForm(const SpacesPtr& spaces,
                     LinearTerms      terms)
    -> LinearForm<SpacesPtr, LinearTerms>
{
  return LinearForm<SpacesPtr, LinearTerms>
                    (spaces, terms);
}

template<class TestSpacesPtr, class LinearTerms, class NewTestSpacesPtr>
auto replaceTestSpaces(
    const LinearForm<TestSpacesPtr, LinearTerms>& linearForm,
    const NewTestSpacesPtr& newTestSpaces)
  -> LinearForm<NewTestSpacesPtr, LinearTerms> {
  return LinearForm<NewTestSpacesPtr, LinearTerms>
                     (newTestSpaces,
                      linearForm.getTerms());
}

} // end namespace Dune

#endif // DUNE_DPG_LINEARFORM_HH
