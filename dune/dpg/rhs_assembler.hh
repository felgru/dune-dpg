// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RHS_ASSEMBLER_HH
#define DUNE_DPG_RHS_ASSEMBLER_HH

#include <tuple>
#include <functional>
#include <memory>
#include <type_traits>

#include "system_assembler.hh"

namespace Dune {

/**
 * \brief This constructs the right hand side vector of a DPG system.
 *
 * \tparam TestSpaces     tuple of test spaces
 */
template<class TestSpaces>
class RhsAssembler
{
public:
  RhsAssembler () = delete;
  /**
   * \brief constructor for RhsAssembler
   *
   * \note For your convenience, use make_RhsAssembler() instead.
   */
  constexpr RhsAssembler (TestSpaces testSpaces)
             : testSpaces(testSpaces)
  { }

  /**
   * \brief Assemble the rhs vector of a DPG system for a given rhs function.
   *
   * Given a tuple of right hand side functions \p volumeTerms,
   * this assembles the right hand side vector of the corresponding
   * DPG system. \p rhs will be overwritten by this function.
   *
   * \param[out] rhs         the rhs vector of the DPG system
   * \param[in]  volumeTerms the rhs functions describing the DPG system
   * \tparam     VolumeTerms a tuple type of rhs functions
   */
  template <class VolumeTerms>
  void assembleRhs(BlockVector<FieldVector<double,1> >& rhs,
                   VolumeTerms&& volumeTerms);

  /**
   * \brief Apply Dirichlet boundary values to the rhs vector
   *
   * \param[in,out] rhs         the rhs vector of the DPG system
   * \param[in] dirichletNodes  true marks the dofs in the Dirichlet boundary
   * \param[in] value           the Dirichlet boundary value
   * \tparam spaceIndex  the index of the test space on which we apply
   *                     the boundary data
   * \tparam ValueType   we take either constants or functions for \p value
   */
  template <size_t spaceIndex, class ValueType>
  void applyDirichletBoundary(BlockVector<FieldVector<double,1> >& rhs,
                              const std::vector<bool>& dirichletNodes,
                              const ValueType& value);

  /**
   * \brief Does exactly what it says on the tin.
   */
  const TestSpaces& getTestSpaces() const
  { return testSpaces; }

private:
  TestSpaces     testSpaces;
};

/**
 * \brief Creates an RhsAssembler,
 *        deducing the target type from the types of arguments.
 *
 * \param testSpaces     a tuple of test spaces
 */
template<class TestSpaces>
RhsAssembler<TestSpaces>
make_RhsAssembler(TestSpaces testSpaces)
{
  return RhsAssembler<TestSpaces>(testSpaces);
}

template<class TestSpaces>
template <class VolumeTerms>
void RhsAssembler<TestSpaces>::
assembleRhs(BlockVector<FieldVector<double,1> >& rhs,
            VolumeTerms&& volumeTerms)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  size_t globalTotalTestSize =
      fold(zip(globalTestSpaceOffsets, testSpaces),
           (size_t)0, globalOffsetHelper());

  rhs.resize(globalTotalTestSize);
  rhs = 0;

  // Views on the FE bases on a single element
  auto testLocalView     = as_vector(transform(testSpaces, getLocalView()));

  auto testLocalIndexSet     = as_vector(transform(testSpaces,
                                                   getLocalIndexSet()));

  for(const auto& e : elements(gridView)) {

    for_each(testLocalView, applyBind<decltype(e)>(e));
    for_each(zip(testLocalIndexSet, testLocalView),
             make_fused_procedure(bindLocalIndexSet()));

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> >
        localRhs[std::tuple_size<
                 typename std::remove_reference<VolumeTerms>::type>::value];

    using RHSZipHelper = vector<decltype(testLocalView)&,
                                decltype(localRhs)&,
                                decltype(volumeTerms)&,
                                decltype(testSpaces)&>;
    for_each(zip_view<RHSZipHelper>(RHSZipHelper(testLocalView,
                                                 localRhs,
                                                 volumeTerms,
                                                 testSpaces)),
             getVolumeTermHelper());

    auto cpr = fused_procedure<localToGlobalRHSCopier<
                    typename std::remove_reference<decltype(rhs)>::type> >
                 (localToGlobalRHSCopier<
                    typename std::remove_reference<decltype(rhs)>::type>(rhs));
    for_each(zip(localRhs,
                 testLocalIndexSet,
                 globalTestSpaceOffsets),
             std::ref(cpr));

  }
}


template<class TestSpaces>
template <size_t spaceIndex, class ValueType>
void RhsAssembler<TestSpaces>::
applyDirichletBoundary(BlockVector<FieldVector<double,1> >& rhs,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& value)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  static_assert(std::is_arithmetic<ValueType>::value,
                "applyDirichletBoundary not implemented for non arithmetic "
                "boundary data types.");

  const size_t spaceSize =
        std::get<spaceIndex>(testSpaces).size();

  size_t globalOffset;
  {
    /* set up global offsets */
    size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
    fold(zip(globalTestSpaceOffsets, testSpaces),
         (size_t)0, globalOffsetHelper());

    globalOffset = globalTestSpaceOffsets[spaceIndex];
  }

  // Set Dirichlet values
  for (size_t i=0; i<spaceSize; i++)
  {
    if (dirichletNodes[i])
    {
      /* TODO: Needs adaptation when value is a function. */
      rhs[globalOffset+i] = value;
    }
  }
}

} // end namespace Dune

#endif // DUNE_DPG_RHS_ASSEMBLER_HH
