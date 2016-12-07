// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_LINEARFORM_HH
#define DUNE_DPG_LINEARFORM_HH

#include <tuple>
#include <vector>
#include <functional>
#include <memory>
#include <type_traits>
#include <cassert>

#include <dune/common/hybridutilities.hh>

#include "assemble_helper.hh"
#include "assemble_types.hh"
#include "linearintegralterm.hh"
#include "type_traits.hh"

namespace Dune {

  /**
   * \brief This class describes a linear form.
   *
   * \tparam Sps             tuple of (test) spaces
   * \tparam LinearTerms     tuple of LinearIntegralTerm
   */
  template<class Sps, class LinearTerms>
  class LinearForm
  {
  public:
    typedef std::decay_t<Sps> Spaces;
    //! tuple type for the local views of the test spaces
    typedef detail::getLocalViews_t<Spaces>  LocalViews;

    LinearForm () = delete;
    /**
     * \brief constructor for LinearForm
     *
     * \note For your convenience, use make_LinearForm() instead.
     */
    constexpr LinearForm (const Spaces&       spaces,
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
    void bind(const LocalViews& lv)
    {
      constexpr bool usesOptimalTestBasis =
            is_OptimalTestSpace<
                typename std::tuple_element<std::tuple_size<Spaces>::value-1,
                                            Spaces>::type
            >::value;

      localViews = std::addressof(lv);

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
    const Spaces& getSpaces() const
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
    Spaces       spaces;
    LinearTerms  terms;

    size_t localSpaceOffsets[std::tuple_size<Spaces>::value];
    size_t localTotalSpaceSize;

    const LocalViews* localViews;
  };

/**
 * \brief Creates a LinearForm for a saddlepoint formulation,
 *        deducing the target type from the types of arguments.
 *
 * \param testSpaces      a tuple of test spaces
 * \param solutionSpaces  a tuple of solution spaces
 * \param terms           a tuple of LinearIntegralTerm
 */
template<class TestSpaces, class SolutionSpaces, class LinearTerms>
auto make_Saddlepoint_LinearForm(TestSpaces      testSpaces,
                                 SolutionSpaces  solutionSpaces,
                                 LinearTerms     terms)
    -> LinearForm<typename std::remove_reference<decltype(
            std::tuple_cat(std::declval<TestSpaces>(),
                           std::declval<SolutionSpaces>())
            )>::type, LinearTerms>
{
  // TODO: This does not really make sense because
  //       of the spaceindex in LinearTerms.
  static_assert(std::tuple_size<SolutionSpaces>::value > 0,
      "Use make_LinearForm and concatenate spaces");
  using Spaces
    = typename std::remove_reference<decltype(
            std::tuple_cat(std::declval<TestSpaces>(),
                           std::declval<SolutionSpaces>()
            ))>::type;
  return LinearForm<Spaces, LinearTerms>
                    (std::tuple_cat(testSpaces, solutionSpaces), terms);
}

/**
 * \brief Creates a LinearForm (i.e. for a saddlepoint formulation),
 *        deducing the target type from the types of arguments.
 *
 * \param spaces          a tuple of spaces
 * \param terms           a tuple of LinearIntegralTerm
 */
template<class Spaces, class LinearTerms>
auto make_LinearForm(Spaces      spaces,
                     LinearTerms     terms)
    -> LinearForm<Spaces, LinearTerms>
{
  return LinearForm<Spaces, LinearTerms>
                    (spaces, terms);
}

template<class TestSpaces, class LinearTerms, class NewTestSpaces>
auto replaceTestSpaces(
    const LinearForm<TestSpaces, LinearTerms>& linearForm,
    NewTestSpaces&& newTestSpaces)
  -> LinearForm<NewTestSpaces, LinearTerms> {
  return LinearForm<NewTestSpaces, LinearTerms>
                     (std::forward<NewTestSpaces>(newTestSpaces),
                      linearForm.getTerms());
}

} // end namespace Dune

#endif // DUNE_DPG_LINEARFORM_HH
