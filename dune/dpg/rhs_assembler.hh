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
 * class RhsAssembler
 *
 * \tparam TestSpaces     tuple of test spaces
 */
template<class TestSpaces>
class RhsAssembler
{
public:
  RhsAssembler () = delete;
  constexpr RhsAssembler (TestSpaces testSpaces)
             : testSpaces(testSpaces)
  { };

  template <class VolumeTerms>
  void assembleRhs(BlockVector<FieldVector<double,1> >& rhs,
                   VolumeTerms&& volumeTerms);

  template <size_t spaceIndex, class ValueType>
  void applyDirichletBoundary(BlockVector<FieldVector<double,1> >& rhs,
                              const std::vector<bool>& dirichletNodes,
                              const ValueType& value);

  const TestSpaces& getTestSpaces() const
  { return testSpaces; };

private:
  TestSpaces     testSpaces;
};

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

  // Get the grid view from the finite element basis
  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  /* TODO: make this a vector of pointers to localVolumeTerm */
  auto localVolumeTerms =
      as_vector(transform(volumeTerms,
                          getLocalVolumeTerm<GridView>(gridView)));

  auto testBasisIndexSet = as_vector(transform(testSpaces, getIndexSet()));

  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  size_t globalTotalTestSize = 0;

  fold(zip(globalTestSpaceOffsets, testBasisIndexSet),
       0, globalOffsetHelper());
  globalTotalTestSize =
      globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value-1]
      + at_c<std::tuple_size<TestSpaces>::value-1>(testBasisIndexSet).size();

  // set rhs to correct length -- the total number of basis vectors in the basis
  rhs.resize(globalTotalTestSize);

  // Set all entries to zero
  rhs = 0;

  // Views on the FE bases on a single element
  auto testLocalView     = as_vector(transform(testSpaces, getLocalView()));

  auto testLocalIndexSet     = as_vector(transform(testBasisIndexSet,
                                                   getLocalIndexSet()));

  // A loop over all elements of the grid
  for(const auto& e : elements(gridView)) {

    // Bind the local FE basis view to the current element
    for_each(testLocalView, applyBind<decltype(e)>(e));

    for_each(zip(testLocalIndexSet, testLocalView),
             make_fused_procedure(bindLocalIndexSet()));

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> >
        localRhs[std::tuple_size<
                 typename std::remove_reference<VolumeTerms>::type>::value];
    for_each(localVolumeTerms, applyBind<decltype(e)>(e));

    using RHSZipHelper = vector<decltype(testLocalView)&,
                                decltype(localRhs)&,
                                decltype(localVolumeTerms)&>;
    for_each(zip_view<RHSZipHelper>(RHSZipHelper(testLocalView,
                                                 localRhs,
                                                 localVolumeTerms)),
             getVolumeTermHelper());

    auto cpr = fused_procedure<localToGlobalRHSCopier<
                        typename remove_reference<decltype(rhs)>::type> >
                    (localToGlobalRHSCopier<
                        typename remove_reference<decltype(rhs)>::type>(rhs));
    for_each(zip(localRhs,
                 testLocalIndexSet,
                 globalTestSpaceOffsets),
             std::ref(cpr));

  }

  /* free memory handled by raw pointers */
  for_each(testLocalIndexSet,     default_deleter());
  for_each(testLocalView,         default_deleter());
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
        std::get<spaceIndex>(testSpaces).indexSet().size();

  size_t globalOffset;
  {
    // Total number of degrees of freedom
    auto testBasisIndexSet = as_vector(transform(testSpaces, getIndexSet()));

    /* set up global offsets */
    size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
    size_t globalTotalTestSize = 0;

    fold(zip(globalTestSpaceOffsets, testBasisIndexSet),
         0, globalOffsetHelper());
    globalTotalTestSize =
        globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value-1]
        + at_c<std::tuple_size<TestSpaces>::value-1>(testBasisIndexSet).size();

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
};

} // end namespace Dune

#endif // DUNE_DPG_RHS_ASSEMBLER_HH
