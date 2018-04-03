// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_BILINEARFORM_HH
#define DUNE_DPG_BILINEARFORM_HH

#include <tuple>
#include <vector>
#include <memory>
#include <type_traits>

#include <boost/hana.hpp>

#include <dune/common/hybridutilities.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include "assemble_helper.hh"
#include "assemble_types.hh"
#include "integralterm.hh"
#include "spacetuple.hh"

namespace Dune {

  /**
   * \brief This class describes a bilinear form.
   *
   * \tparam TSpaces         shared_ptr of tuple of test spaces
   * \tparam SolSpaces       shared_ptr of tuple of solution spaces
   * \tparam BilinearTerms   tuple of IntegralTerm
   */
  template<class TSpaces, class SolSpaces, class BilinearTerms>
  class BilinearForm
  {
    static_assert(is_SpaceTuplePtr<TSpaces>::value,
        "TSpaces needs to be a SpaceTuplePtr!");
    static_assert(is_SpaceTuplePtr<SolSpaces>::value,
        "SolSpaces needs to be a SpaceTuplePtr!");
  public:
    //! shared pointer to tuple type of test spaces
    using TestSpacesPtr = TSpaces;
    //! shared pointer to tuple type of solution spaces
    using SolutionSpacesPtr = SolSpaces;
    //! tuple type of test spaces
    using TestSpaces = typename TestSpacesPtr::element_type;
    //! tuple type of solution spaces
    using SolutionSpaces = typename SolutionSpacesPtr::element_type;
    //! tuple type of bilinear form terms
    using Terms = BilinearTerms;
    //! tuple type for the local views of the test spaces
    using TestLocalViews = detail::getLocalViews_t<TestSpaces>;
    //! tuple type for the local views of the solution spaces
    using SolutionLocalViews = detail::getLocalViews_t<SolutionSpaces>;

    BilinearForm () = delete;
    /**
     * \brief constructor for BilinearForm
     *
     * \note For your convenience, use make_BilinearForm() instead.
     */
    constexpr BilinearForm (const TestSpacesPtr&     testSpaces,
                            const SolutionSpacesPtr& solutionSpaces,
                            const BilinearTerms&  terms)
               : testSpaces(testSpaces),
                 solutionSpaces(solutionSpaces),
                 terms(terms),
                 testLocalViews(nullptr),
                 solutionLocalViews(nullptr)
    { }

    /**
     * \brief Compute the stiffness matrix for a single element.
     *
     * \pre bind() has to be called on the current element before.
     *
     * \param[out] elementMatrix will hold the local system matrix
     */
    template <class MatrixType>
    void getLocalMatrix(MatrixType& elementMatrix) const
    {
      // Set all matrix entries to zero
      elementMatrix.setSize(localTotalTestSize,
                            localTotalSolutionSize);
      elementMatrix = 0;

      Hybrid::forEach(terms,
              detail::getLocalMatrixHelper
                      <MatrixType,
                       TestSpaces,
                       SolutionSpaces>
                      (*testLocalViews,
                       *solutionLocalViews,
                       elementMatrix,
                       localTestSpaceOffsets,
                       localSolutionSpaceOffsets));
    }

    /**
     * \brief Creates the occupation pattern for the system matrix.
     *
     * This creates an index set of all those entries of the system
     * matrix that might be non-zero.
     *
     * \param nb  the matrix index set of the non-zero entries
     * \param testShift      offset for the test space indices
     * \param solutionShift  offset for the solution space indices
     *
     * \tparam mirror  if true, the indices additionally get written mirrored,
     *                 this is used for the saddlepoint formulation
     */
    template<bool mirror=false>
    void getOccupationPattern(MatrixIndexSet& nb,
                              size_t testShift, size_t solutionShift) const;

    /**
     * \brief Bind the BilinearForm to the given localViews.
     *
     * \pre The given localViews have to be bound to the same element.
     */
    void bind(const TestLocalViews& tlv, const SolutionLocalViews& slv)
    {
      testLocalViews     = std::addressof(tlv);
      solutionLocalViews = std::addressof(slv);

      using namespace Dune::detail;

      /* set up local offsets */
      localTotalTestSize = computeOffsets(localTestSpaceOffsets,
                                          tlv);

      localTotalSolutionSize = computeOffsets(localSolutionSpaceOffsets,
                                              slv);
    }

    /**
     * \brief Does exactly what it says on the tin.
     */
    TestSpacesPtr getTestSpaces() const
    { return testSpaces; }

    /**
     * \brief Does exactly what it says on the tin.
     */
    SolutionSpacesPtr getSolutionSpaces() const
    { return solutionSpaces; }

    /**
     * \brief Does exactly what it says on the tin.
     */
    const BilinearTerms& getTerms() const
    {
      return terms;
    }

    using TestSpaceIndexArray = size_t[std::tuple_size<TestSpaces>::value];
    using SolutionSpaceIndexArray
        = size_t[std::tuple_size<SolutionSpaces>::value];

    /**
     * \brief Does exactly what it says on the tin.
     */
    const TestSpaceIndexArray& getLocalTestSpaceOffsets() const
    { return localTestSpaceOffsets; }

    /**
     * \brief Does exactly what it says on the tin.
     */
    const SolutionSpaceIndexArray& getLocalSolutionSpaceOffsets() const
    { return localSolutionSpaceOffsets; }

  private:
    TestSpacesPtr     testSpaces;
    SolutionSpacesPtr solutionSpaces;
    BilinearTerms     terms;

    size_t localTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
    size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
    size_t localTotalTestSize,
           localTotalSolutionSize;

    const TestLocalViews*     testLocalViews;
    const SolutionLocalViews* solutionLocalViews;
  };

/**
 * \brief Creates a BilinearForm,
 *        deducing the target type from the types of arguments.
 *
 * \param testSpaces     a shared_ptr to a tuple of test spaces
 * \param solutionSpaces a shared_ptr to a tuple of solution spaces
 * \param terms          a tuple of IntegralTerm
 */
template<class TestSpacesPtr, class SolutionSpacesPtr, class BilinearTerms>
auto make_BilinearForm(const TestSpacesPtr&     testSpaces,
                       const SolutionSpacesPtr& solutionSpaces,
                       BilinearTerms            terms)
    -> BilinearForm<TestSpacesPtr, SolutionSpacesPtr, BilinearTerms>
{
  return BilinearForm<TestSpacesPtr, SolutionSpacesPtr, BilinearTerms>
                      (testSpaces,
                       solutionSpaces,
                       terms);
}

template<class TestSpacesPtr, class SolutionSpacesPtr, class BilinearTerms,
         class NewTestSpacesPtr>
auto replaceTestSpaces(
    const BilinearForm<TestSpacesPtr, SolutionSpacesPtr, BilinearTerms>&
                                                              bilinearForm,
    const NewTestSpacesPtr& newTestSpaces) {
  return make_BilinearForm
                     (newTestSpaces,
                      bilinearForm.getSolutionSpaces(),
                      bilinearForm.getTerms());
}


template<class TestSpacesPtr, class SolutionSpacesPtr, class BilinearTerms>
template<bool mirror>
void BilinearForm<TestSpacesPtr, SolutionSpacesPtr, BilinearTerms>::
getOccupationPattern(MatrixIndexSet& nb,
                     const size_t testShift,
                     const size_t solutionShift) const
{
  using namespace Dune::detail;

  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  computeOffsets(globalTestSpaceOffsets, *testSpaces, testShift);
  computeOffsets(globalSolutionSpaceOffsets, *solutionSpaces, solutionShift);

  // A view on the FE basis on a single element
  auto solutionLocalViews = getLocalViews(*solutionSpaces);
  auto testLocalViews     = getLocalViews(*testSpaces);

  auto solutionLocalIndexSets = getLocalIndexSets(*solutionSpaces);
  auto testLocalIndexSets = getLocalIndexSets(*testSpaces);

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  const GridView gridView = std::get<0>(*testSpaces).gridView();

  /* create set of index pairs from bilinearTerms to loop over. */
  namespace hana = boost::hana;
  auto indexPairs = hana::to<hana::set_tag>(
      hana::transform(hana::to<hana::tuple_tag>(
          hana::make_range(hana::int_c<0>,
            hana::int_c<std::tuple_size<BilinearTerms>::value>)),
        [](auto i) {
          using Term = std::tuple_element_t<i.value, BilinearTerms>;
          return hana::tuple<std::tuple_element_t<0, Term>,
                             std::tuple_element_t<1, Term>>{};
        }));
  using IndexPairs = decltype(hana::to<hana::tuple_tag>(indexPairs));

  for(const auto& e : elements(gridView))
  {
    bindLocalViews(solutionLocalViews, e);
    bindLocalViews(testLocalViews, e);

    bindLocalIndexSets(solutionLocalIndexSets, solutionLocalViews);
    bindLocalIndexSets(testLocalIndexSets, testLocalViews);

    detail::getOccupationPattern<IndexPairs, mirror>
                        (testLocalIndexSets,
                         solutionLocalIndexSets,
                         globalTestSpaceOffsets,
                         globalSolutionSpaceOffsets,
                         nb);
  }
}

} // end namespace Dune

#endif // DUNE_DPG_BILINEARFORM_HH
