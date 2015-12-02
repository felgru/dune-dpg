// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_LINEARFORM_HH
#define DUNE_DPG_LINEARFORM_HH

#include <tuple>
#include <vector>
#include <functional>
#include <type_traits>
#include <cassert>

#include <dune/dpg/common/memory.hh>

#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/adapted/array.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/container/set/convert.hpp>
#include <boost/fusion/sequence/intrinsic/size.hpp>
#include <boost/fusion/algorithm/auxiliary/copy.hpp>
#include <boost/fusion/algorithm/transformation/transform.hpp>
#include <boost/fusion/algorithm/transformation/zip.hpp>
#include <boost/fusion/algorithm/iteration/accumulate.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/functional/generation/make_fused_procedure.hpp>

#include "assemble_helper.hh"
#include "assemble_types.hh"
#include "linearintegralterm.hh"

namespace Dune {

  /**
   * \brief This class describes a linear form.
   *
   * \tparam Sps             tuple of (test) spaces
   * \tparam LinearTerms     tuple of LinearIntegralTerm
   * \tparam FormulationType either SaddlepointFormulation or DPGFormulation
   */
  template<class Sps, class LinearTerms, class FormulationType>
  class LinearForm
  {
  public:
    typedef std::decay_t<Sps> Spaces;
    //! tuple type for the local views of the test spaces
    typedef typename boost::fusion::result_of::as_vector<
        typename boost::fusion::result_of::
          transform<Spaces, detail::getLocalView>::type
        >::type LocalViews;

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
      elementVector.resize(localTotalSize);
      elementVector = 0;

      boost::fusion::for_each(terms,
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
      constexpr bool isSaddlepoint =
                std::is_same<
                typename std::decay<FormulationType>::type
              , SaddlepointFormulation
            >::value;
      //static_assert(!isSaddlepoint,
      //              "LinearForm only implemented for DPGFormulation.");

      localViews = std::addressof(lv);

      using namespace boost::fusion;
      using namespace Dune::detail;

      constexpr size_t size = result_of::size<LocalViews>::value;

      /* set up local offsets */
      if(isSaddlepoint) {
        fold(zip(localSpaceOffsets, lv), (size_t)0, offsetHelper());
      } else { /* DPG formulation */
        for(size_t i=0; i<size; ++i)
        {
          localSpaceOffsets[i] = 0;
        }
      }
      localTotalSize =
          localSpaceOffsets[size-1]
          + at_c<size-1>(lv).size();
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
    size_t localTotalSize;

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
            )>::type, LinearTerms, SaddlepointFormulation>
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
  return LinearForm<Spaces, LinearTerms, SaddlepointFormulation>
                    (std::tuple_cat(testSpaces, solutionSpaces), terms);
}

/**
 * \brief Creates a LinearForm for a DPG formulation,
 *        deducing the target type from the types of arguments.
 *
 * \param testSpaces  a tuple of test spaces
 * \param terms       a tuple of LinearIntegralTerm
 */
template<class TestSpaces, class LinearTerms>
auto make_DPG_LinearForm(TestSpaces   testSpaces,
                         LinearTerms  terms)
    -> LinearForm<TestSpaces, LinearTerms, DPGFormulation>
{
  return LinearForm<TestSpaces, LinearTerms, DPGFormulation>
                    (testSpaces, terms);
}

/**
 * \brief Creates a LinearForm for a DPG formulation consistent with DPGSystemAssembler,
 *        deducing the target type from the types of arguments.
 *
 * \param testSpaces  a tuple of test spaces
 * \param terms       a tuple of LinearIntegralTerm
 */
template<class TestSpaces, class LinearTerms>
auto make_DPGLinearForm(TestSpaces   testSpaces,
                        LinearTerms  terms)
    -> LinearForm<TestSpaces, LinearTerms, SaddlepointFormulation>
{
  return LinearForm<TestSpaces, LinearTerms, SaddlepointFormulation>
                    (testSpaces, terms);
}

/**
 * \brief Creates a LinearForm (i.e. for a saddlepoint formulation),
 *        deducing the target type from the types of arguments.
 *
 * \param Spaces          a tuple of spaces
 * \param terms           a tuple of LinearIntegralTerm
 */
template<class Spaces, class LinearTerms>
auto make_LinearForm(Spaces      spaces,
                     LinearTerms     terms)
    -> LinearForm<Spaces, LinearTerms, SaddlepointFormulation>
{
  return LinearForm<Spaces, LinearTerms, SaddlepointFormulation>
                    (spaces, terms);
}

template<class TestSpaces, class LinearTerms,
         class FormulationType, class NewTestSpaces>
auto replaceTestSpaces(
    const LinearForm<TestSpaces, LinearTerms, FormulationType>& linearForm,
    NewTestSpaces&& newTestSpaces)
  -> LinearForm<NewTestSpaces, LinearTerms, FormulationType> {
  return LinearForm<NewTestSpaces, LinearTerms, FormulationType>
                     (std::forward<NewTestSpaces>(newTestSpaces),
                      linearForm.getTerms());
}

} // end namespace Dune

#endif // DUNE_DPG_LINEARFORM_HH
