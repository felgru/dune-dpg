// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_INNERPRODUCT_HH
#define DUNE_DPG_INNERPRODUCT_HH

#include <tuple>
#include <vector>
#include <functional>
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

namespace Dune {

  /**
   * \brief This class describes an inner product.
   *
   * \tparam TSpaces            tuple of test spaces
   * \tparam InnerProductTerms  tuple of IntegralTerm
   */
  template<class TSpaces, class InnerProductTerms>
  class InnerProduct
  {
  public:
    //! tuple type of test spaces
    typedef TSpaces TestSpaces;
    //! tuple type of inner product terms
    typedef InnerProductTerms Terms;
    //! tuple type for the local views of the test spaces
    using TestLocalViews = detail::getLocalViews_t<TestSpaces>;

    InnerProduct () = delete;
    /**
     * \brief constructor for InnerProduct
     *
     * \note For your convenience, use make_InnerProduct() instead.
     */
    constexpr InnerProduct (const std::shared_ptr<TestSpaces>& testSpaces,
                            const InnerProductTerms&           terms)
               : testSpaces(testSpaces),
                 terms(terms),
                 testLocalViews(nullptr)
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
                            localTotalTestSize);
      elementMatrix = 0;

      Hybrid::forEach(terms,
              detail::getLocalMatrixHelper
                      <MatrixType,
                       TestSpaces,
                       TestSpaces>
                      (*testLocalViews,
                       *testLocalViews,
                       elementMatrix,
                       localTestSpaceOffsets,
                       localTestSpaceOffsets));
    }

    /**
     * \brief Creates the occupation pattern for the system matrix.
     *
     * This creates an index set of all those entries of the system
     * matrix that might be non-zero.
     *
     * \param nb  the matrix index set of the non-zero entries
     */
    void getOccupationPattern(MatrixIndexSet& nb) const;

    /**
     * \brief Bind the InnerProduct to the given localViews.
     *
     * \pre The given localViews have to be bound to the same element.
     */
    void bind(const TestLocalViews& tlv)
    {
      testLocalViews = std::addressof(tlv);

      /* set up local offsets */
      localTotalTestSize = detail::computeOffsets(localTestSpaceOffsets, tlv);
    }

    /**
     * \brief Does exactly what it says on the tin.
     */
    std::shared_ptr<TestSpaces> getTestSpaces() const
    { return testSpaces; }

    /**
     * \brief Does exactly what it says on the tin.
     */
    const InnerProductTerms& getTerms() const
    {
      return terms;
    }

  private:
    std::shared_ptr<TestSpaces> testSpaces;
    InnerProductTerms           terms;
    size_t localTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
    size_t localTotalTestSize;
    const TestLocalViews* testLocalViews;
  };

/**
 * \brief Creates an InnerProduct,
 *        deducing the target type from the types of arguments.
 *
 * \param testSpaces  a shared_ptr to a tuple of test spaces
 * \param terms       a tuple of IntegralTerm
 */
template<class TestSpaces, class InnerProductTerms>
auto make_InnerProduct(const std::shared_ptr<TestSpaces>& testSpaces,
                       InnerProductTerms                  terms)
    -> InnerProduct<TestSpaces, InnerProductTerms>
{
  return InnerProduct<TestSpaces, InnerProductTerms> (testSpaces, terms);
}

template<class TestSpaces, class InnerProductTerms, class NewTestSpaces>
auto replaceTestSpaces(
    const InnerProduct<TestSpaces, InnerProductTerms>& innerProduct,
    const shared_ptr<NewTestSpaces>& newTestSpaces) {
  return make_InnerProduct(newTestSpaces,
                           innerProduct.getTerms());
}

template<class TestSpaces, class InnerProductTerms>
void InnerProduct<TestSpaces, InnerProductTerms>::
getOccupationPattern(MatrixIndexSet& nb) const
{
  using namespace Dune::detail;

  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  computeOffsets(globalTestSpaceOffsets, *testSpaces);

  auto testLocalViews     = getLocalViews(*testSpaces);
  auto testLocalIndexSets = getLocalIndexSets(*testSpaces);

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(*testSpaces).gridView();

  /* create set of index pairs from innerProductTerms to loop over. */
  namespace hana = boost::hana;
  auto indexPairs = hana::to<hana::set_tag>(
      hana::transform(hana::to<hana::tuple_tag>(
          hana::make_range(hana::int_c<0>,
            hana::int_c<std::tuple_size<InnerProductTerms>::value>)),
        [](auto i) {
          using Term = std::tuple_element_t<i.value, InnerProductTerms>;
          return hana::tuple<std::tuple_element_t<0, Term>,
                             std::tuple_element_t<1, Term>>{};
        }));
  using IndexPairs = decltype(hana::to<hana::tuple_tag>(indexPairs));

  for(const auto& e : elements(gridView))
  {
    bindLocalViews(testLocalViews, e);
    bindLocalIndexSets(testLocalIndexSets, testLocalViews);

    detail::getOccupationPattern<IndexPairs, false>
                        (testLocalIndexSets,
                         testLocalIndexSets,
                         globalTestSpaceOffsets,
                         globalTestSpaceOffsets,
                         nb);
  }
}

} // end namespace Dune

#endif // DUNE_DPG_INNERPRODUCT_HH
