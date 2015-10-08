// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_BILINEARFORM_HH
#define DUNE_DPG_BILINEARFORM_HH

#include <tuple>
#include <vector>
#include <functional>
#include <memory>
#include <type_traits>

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

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include "assemble_helper.hh"
#include "assemble_types.hh"
#include "integralterm.hh"

namespace Dune {

  /**
   * \brief This class describes a bilinear form.
   *
   * \tparam TSpaces         tuple of test spaces
   * \tparam SolSpaces       tuple of solution spaces
   * \tparam BilinearTerms   tuple of IntegralTerm
   * \tparam FormulationType either SaddlepointFormulation or DPGFormulation
   */
  template<class TSpaces, class SolSpaces, class BilinearTerms, class FormulationType>
  class BilinearForm
  {
  public:
    //! tuple type of test spaces
    typedef TSpaces TestSpaces;
    //! tuple type of solution spaces
    typedef SolSpaces SolutionSpaces;
    //! tuple type for the local views of the test spaces
    typedef typename boost::fusion::result_of::as_vector<
        typename boost::fusion::result_of::
        transform<TestSpaces, detail::getLocalView>::type>::type TestLocalView;
    //! tuple type for the local views of the solution spaces
    typedef typename boost::fusion::result_of::as_vector<
        typename boost::fusion::result_of::
        transform<SolutionSpaces, detail::getLocalView>::type
        >::type SolutionLocalView;

    BilinearForm () = delete;
    /**
     * \brief constructor for BilinearForm
     *
     * \note For your convenience, use make_BilinearForm() instead.
     */
    constexpr BilinearForm (const TestSpaces&     testSpaces,
                            const SolutionSpaces& solutionSpaces,
                            const BilinearTerms&  terms)
               : testSpaces(testSpaces),
                 solutionSpaces(solutionSpaces),
                 terms(terms),
                 testLocalView(nullptr),
                 solutionLocalView(nullptr)
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

      boost::fusion::for_each(terms,
              detail::getLocalMatrixHelper
                      <MatrixType,
                       TestSpaces,
                       SolutionSpaces>
                      (*testLocalView,
                       *solutionLocalView,
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
    void bind(const TestLocalView& tlv, const SolutionLocalView& slv)
    {
      constexpr bool isSaddlepoint =
                std::is_same<
                typename std::decay<FormulationType>::type
              , SaddlepointFormulation
            >::value;

      testLocalView     = std::addressof(tlv);
      solutionLocalView = std::addressof(slv);

      using namespace boost::fusion;
      using namespace Dune::detail;

      using TestSize     = typename result_of::size<TestLocalView>::type;
      using SolutionSize = typename result_of::size<SolutionLocalView>::type;

      /* set up local offsets */
      if(isSaddlepoint) {
        fold(zip(localTestSpaceOffsets, tlv), (size_t)0, offsetHelper());
      } else { /* DPG formulation */
        for(size_t i=0; i<std::tuple_size<TestSpaces>::value; ++i)
        {
          localTestSpaceOffsets[i] = 0;
        }
      }
      localTotalTestSize =
          localTestSpaceOffsets[TestSize::value-1]
          + at_c<TestSize::value-1>(tlv).size();

      localTotalSolutionSize =
          fold(zip(localSolutionSpaceOffsets, slv),
               (size_t)0, offsetHelper());
    }

    /**
     * \brief Does exactly what it says on the tin.
     */
    const TestSpaces& getTestSpaces() const
    { return testSpaces; }

    /**
     * \brief Does exactly what it says on the tin.
     */
    const SolutionSpaces& getSolutionSpaces() const
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
    TestSpaces     testSpaces;
    SolutionSpaces solutionSpaces;
    BilinearTerms  terms;

    size_t localTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
    size_t localSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
    size_t localTotalTestSize,
           localTotalSolutionSize;

    const TestLocalView*     testLocalView;
    const SolutionLocalView* solutionLocalView;
  };

/**
 * \brief Creates a BilinearForm,
 *        deducing the target type from the types of arguments.
 *
 * \param testSpaces     a tuple of test spaces
 * \param solutionSpaces a tuple of solution spaces
 * \param terms          a tuple of IntegralTerm
 */
template<class TestSpaces, class SolutionSpaces, class BilinearTerms>
auto make_BilinearForm(TestSpaces     testSpaces,
                       SolutionSpaces solutionSpaces,
                       BilinearTerms  terms)
    -> BilinearForm<TestSpaces, SolutionSpaces, BilinearTerms, SaddlepointFormulation>
{
  return BilinearForm<TestSpaces, SolutionSpaces, BilinearTerms, SaddlepointFormulation>
                      (testSpaces,
                       solutionSpaces,
                       terms);
}


namespace detail {

  template<class FormulationType, class TestSpaces, class SolutionSpaces, class BilinearTerms>
  auto make_BilinearForm(TestSpaces     testSpaces,
                         SolutionSpaces solutionSpaces,
                         BilinearTerms  terms)
      -> BilinearForm<TestSpaces, SolutionSpaces, BilinearTerms, FormulationType>
  {
    return BilinearForm<TestSpaces, SolutionSpaces, BilinearTerms, FormulationType>
                        (testSpaces,
                         solutionSpaces,
                         terms);
  }

}


template<class TestSpaces, class SolutionSpaces, class BilinearTerms, class FormulationType>
template<bool mirror>
void BilinearForm<TestSpaces, SolutionSpaces, BilinearTerms, FormulationType>::
getOccupationPattern(MatrixIndexSet& nb, size_t testShift, size_t solutionShift) const
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  constexpr bool isSaddlepoint =
        std::is_same<
             typename std::decay<FormulationType>::type
           , SaddlepointFormulation
        >::value;


  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  if(isSaddlepoint) {
    fold(zip(globalTestSpaceOffsets, testSpaces),
         testShift, globalOffsetHelper());
  } else { /* DPG formulation */
    for(size_t i=0; i<std::tuple_size<TestSpaces>::value; ++i)
    {
      globalTestSpaceOffsets[i] = 0;
    }
  }
  fold(zip(globalSolutionSpaceOffsets, solutionSpaces),
       solutionShift, globalOffsetHelper());

  // A view on the FE basis on a single element
  auto solutionLocalView = as_vector(transform(solutionSpaces,
                                               getLocalView()));
  auto testLocalView     = as_vector(transform(testSpaces,
                                               getLocalView()));

  auto solutionLocalIndexSet = as_vector(transform(solutionSpaces,
                                                   getLocalIndexSet()));
  auto testLocalIndexSet     = as_vector(transform(testSpaces,
                                                   getLocalIndexSet()));

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  /* create set of index pairs from bilinearTerms to loop over. */
  typedef typename boost::mpl::fold<
      typename boost::mpl::transform<
          /* This as_vector is probably not needed for boost::fusion 1.58
           * or higher. */
          typename result_of::as_vector<BilinearTerms>::type
        , mpl::firstTwo<boost::mpl::_1>
        >::type
    , boost::mpl::set0<>
    , boost::mpl::insert<boost::mpl::_1,boost::mpl::_2>
    >::type IndexPairs;
  auto indexPairs = IndexPairs{};

  for(const auto& e : elements(gridView))
  {
    for_each(solutionLocalView, applyBind<decltype(e)>(e));
    for_each(testLocalView, applyBind<decltype(e)>(e));

    for_each(zip(solutionLocalIndexSet, solutionLocalView),
             make_fused_procedure(bindLocalIndexSet()));
    for_each(zip(testLocalIndexSet, testLocalView),
             make_fused_procedure(bindLocalIndexSet()));

    auto gOPH = getOccupationPatternHelper<decltype(testLocalView),
                                           decltype(solutionLocalView),
                                           decltype(testLocalIndexSet),
                                           decltype(solutionLocalIndexSet),
                                           mirror>
                        (testLocalView,
                         solutionLocalView,
                         testLocalIndexSet,
                         solutionLocalIndexSet,
                         globalTestSpaceOffsets,
                         globalSolutionSpaceOffsets,
                         nb);
    for_each(indexPairs,
        std::ref(gOPH));
  }
}

} // end namespace Dune

#endif // DUNE_DPG_BILINEARFORM_HH
