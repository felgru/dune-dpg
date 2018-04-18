// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_FUNCTIONS_LOCALINDEXSETITERATION_HH
#define DUNE_DPG_FUNCTIONS_LOCALINDEXSETITERATION_HH

#include <type_traits>

#include <dune/common/version.hh>

#include <dune/dpg/functions/concepts.hh>
#include <dune/functions/functionspacebases/concepts.hh>

namespace Dune {
namespace Functions {

namespace detail {

  template<class LocalIndexSet,
           class UnconstrainedEvaluation,
           class ConstrainedSetUp,
           class ConstrainedEvaluation,
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
           typename std::enable_if<models<Concept::LocalView<typename
                                std::decay_t<LocalIndexSet>::GlobalBasis>,
                          std::decay_t<LocalIndexSet>>()>::type* = nullptr
#else
           typename std::enable_if<models<Concept::LocalIndexSet<typename
                                std::decay_t<LocalIndexSet>::LocalView>,
                          std::decay_t<LocalIndexSet>>()>::type* = nullptr
#endif
                            >
  void iterateOverLocalIndices_impl(
      LocalIndexSet&& localIndexSet,
      UnconstrainedEvaluation&& unconstrainedEvaluation,
      ConstrainedSetUp&&,
      ConstrainedEvaluation&&) {
    using size_type = typename std::decay_t<LocalIndexSet>::size_type;
    for(size_type i=0, i_max=localIndexSet.size(); i<i_max; i++) {
      unconstrainedEvaluation(i, localIndexSet.index(i));
    }
  }

  template<class LocalIndexSet,
           class UnconstrainedEvaluation,
           class ConstrainedSetUp,
           class ConstrainedEvaluation,
#if DUNE_VERSION_NEWER(DUNE_FUNCTIONS,2,7)
           typename std::enable_if<models<Concept::ConstrainedLocalView<
                          typename std::decay_t<LocalIndexSet>::GlobalBasis>,
                        std::decay_t<LocalIndexSet>>()>::type* = nullptr>
#else
           typename std::enable_if<models<Concept::ConstrainedLocalIndexSet<
                          typename std::decay_t<LocalIndexSet>::LocalView>,
                        std::decay_t<LocalIndexSet>>()>::type* = nullptr>
#endif
  void iterateOverLocalIndices_impl(
      LocalIndexSet&& localIndexSet,
      UnconstrainedEvaluation&& unconstrainedEvaluation,
      ConstrainedSetUp&& constrainedSetUp,
      ConstrainedEvaluation&& constrainedEvaluation) {
    using size_type = typename std::decay_t<LocalIndexSet>::size_type;
    auto globalIndex = localIndexSet.indicesLocalGlobal().begin();
    assert(!localIndexSet.indicesLocalGlobal().empty());
    const size_type numConstraints = localIndexSet.constraintsSize();
    size_type i = 0;
    for(size_type c = 0; c < numConstraints; c++) {
      const size_type nextConstraint = i + localIndexSet.constraintOffset(c);
      for (; i < nextConstraint; ++i) {
        unconstrainedEvaluation(i, *(globalIndex++));
      }
      constrainedSetUp(i);
      for(auto w : localIndexSet.constraintWeights(c)) {
        constrainedEvaluation(i, *(globalIndex++), w);
      }
      ++i;
    }
    for (; i < localIndexSet.size(); ++i)
      unconstrainedEvaluation(i, *(globalIndex++));
    assert(globalIndex == localIndexSet.indicesLocalGlobal().end());
  }
}

template<class LocalIndexSet,
         class UnconstrainedEvaluation,
         class ConstrainedSetUp,
         class ConstrainedEvaluation>
void iterateOverLocalIndices(
    LocalIndexSet&& localIndexSet,
    UnconstrainedEvaluation&& unconstrainedEvaluation,
    ConstrainedSetUp&& constrainedSetUp,
    ConstrainedEvaluation&& constrainedEvaluation) {
  detail::iterateOverLocalIndices_impl(
      std::forward<LocalIndexSet>(localIndexSet),
      std::forward<UnconstrainedEvaluation>(unconstrainedEvaluation),
      std::forward<ConstrainedSetUp>(constrainedSetUp),
      std::forward<ConstrainedEvaluation>(constrainedEvaluation));
}

template<class TestLocalIndexSet,
         class SolutionLocalIndexSet,
         class GetLocalMatrixEntry,
         class GetGlobalMatrixEntry>
void addToGlobalMatrix(
    TestLocalIndexSet&& testLocalIndexSet,
    SolutionLocalIndexSet&& solutionLocalIndexSet,
    GetLocalMatrixEntry&& getLocalMatrixEntry,
    GetGlobalMatrixEntry&& getGlobalMatrixEntry) {
  using TestMultiIndex = typename std::decay_t<TestLocalIndexSet>::MultiIndex;
  using SolutionMultiIndex
      = typename std::decay_t<SolutionLocalIndexSet>::MultiIndex;
  iterateOverLocalIndices(
    std::forward<TestLocalIndexSet>(testLocalIndexSet),
    [&](size_t i, TestMultiIndex gi)
    {
      iterateOverLocalIndices(
        solutionLocalIndexSet,
        [&](size_t j, SolutionMultiIndex gj)
        {
          getGlobalMatrixEntry(gi, gj) += getLocalMatrixEntry(i, j);
        },
        [](size_t /* j */){},
        [&](size_t j, TestMultiIndex gj, double wj)
        {
          getGlobalMatrixEntry(gi, gj) += wj * getLocalMatrixEntry(i, j);
        }
      );
    },
    [](size_t /* i */){},
    [&](size_t i, TestMultiIndex gi, double wi)
    {
      iterateOverLocalIndices(
        solutionLocalIndexSet,
        [&](size_t j, SolutionMultiIndex gj)
        {
          getGlobalMatrixEntry(gi, gj) += wi * getLocalMatrixEntry(i, j);
        },
        [](size_t /* j */){},
        [&](size_t j, TestMultiIndex gj, double wj)
        {
          getGlobalMatrixEntry(gi, gj) += wi * wj * getLocalMatrixEntry(i, j);
        }
      );
    }
  );
}

} // namespace Dune::Functions
} // namespace Dune

#endif // DUNE_DPG_FUNCTIONS_LOCALINDEXSETITERATION_HH
