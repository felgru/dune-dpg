// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_FUNCTIONS_LOCALINDEXSETITERATION_HH
#define DUNE_DPG_FUNCTIONS_LOCALINDEXSETITERATION_HH

#include <type_traits>

#include <dune/dpg/functions/concepts.hh>
#include <dune/functions/functionspacebases/concepts.hh>

namespace Dune {
namespace Functions {

namespace detail {

  template<class LocalView,
           class UnconstrainedEvaluation,
           class ConstrainedSetUp,
           class ConstrainedEvaluation,
           typename std::enable_if<models<Concept::LocalView<typename
                                std::decay_t<LocalView>::GlobalBasis>,
                          std::decay_t<LocalView>>()>::type* = nullptr>
  void iterateOverLocalIndices_impl(
      LocalView&& localView,
      UnconstrainedEvaluation&& unconstrainedEvaluation,
      ConstrainedSetUp&&,
      ConstrainedEvaluation&&) {
    using size_type = typename std::decay_t<LocalView>::size_type;
    for(size_type i=0, i_max=localView.size(); i<i_max; i++) {
      unconstrainedEvaluation(i, localView.index(i));
    }
  }

  template<class LocalView,
           class UnconstrainedEvaluation,
           class ConstrainedSetUp,
           class ConstrainedEvaluation,
           typename std::enable_if<models<Concept::ConstrainedLocalView<
                          typename std::decay_t<LocalView>::GlobalBasis>,
                        std::decay_t<LocalView>>()>::type* = nullptr>
  void iterateOverLocalIndices_impl(
      LocalView&& localView,
      UnconstrainedEvaluation&& unconstrainedEvaluation,
      ConstrainedSetUp&& constrainedSetUp,
      ConstrainedEvaluation&& constrainedEvaluation) {
    using size_type = typename std::decay_t<LocalView>::size_type;
    auto globalIndex = localView.indicesLocalGlobal().begin();
    assert(!localView.indicesLocalGlobal().empty());
    const size_type numConstraints = localView.constraintsSize();
    size_type i = 0;
    for(size_type c = 0; c < numConstraints; c++) {
      const size_type nextConstraint = i + localView.constraintOffset(c);
      for (; i < nextConstraint; ++i) {
        unconstrainedEvaluation(i, *(globalIndex++));
      }
      constrainedSetUp(i);
      for(auto w : localView.constraintWeights(c)) {
        constrainedEvaluation(i, *(globalIndex++), w);
      }
      ++i;
    }
    for (; i < localView.size(); ++i)
      unconstrainedEvaluation(i, *(globalIndex++));
    assert(globalIndex == localView.indicesLocalGlobal().end());
  }
}

template<class LocalView,
         class UnconstrainedEvaluation,
         class ConstrainedSetUp,
         class ConstrainedEvaluation>
void iterateOverLocalIndices(
    LocalView&& localView,
    UnconstrainedEvaluation&& unconstrainedEvaluation,
    ConstrainedSetUp&& constrainedSetUp,
    ConstrainedEvaluation&& constrainedEvaluation) {
  detail::iterateOverLocalIndices_impl(
      std::forward<LocalView>(localView),
      std::forward<UnconstrainedEvaluation>(unconstrainedEvaluation),
      std::forward<ConstrainedSetUp>(constrainedSetUp),
      std::forward<ConstrainedEvaluation>(constrainedEvaluation));
}

template<class LocalView,
         class GlobalVector,
         class LocalVector>
void copyToLocalVector(const GlobalVector& globalVector,
                       LocalVector& localVector,
                       LocalView&& localView) {
  localVector.resize(localView.size());

  iterateOverLocalIndices(
    localView,
    [&](size_t j, auto gj) {
      localVector[j] = globalVector[gj[0]];
    },
    [&](size_t j) { localVector[j] = 0; },
    [&](size_t j, auto gj, double wj) {
      localVector[j] += wj * globalVector[gj[0]];
    }
  );
}

template<class TestLocalView,
         class SolutionLocalView,
         class GetLocalMatrixEntry,
         class GetGlobalMatrixEntry>
void addToGlobalMatrix(
    TestLocalView&& testLocalView,
    SolutionLocalView&& solutionLocalView,
    GetLocalMatrixEntry&& getLocalMatrixEntry,
    GetGlobalMatrixEntry&& getGlobalMatrixEntry) {
  using TestMultiIndex = typename std::decay_t<TestLocalView>::MultiIndex;
  using SolutionMultiIndex
      = typename std::decay_t<SolutionLocalView>::MultiIndex;
  iterateOverLocalIndices(
    std::forward<TestLocalView>(testLocalView),
    [&](size_t i, TestMultiIndex gi)
    {
      iterateOverLocalIndices(
        solutionLocalView,
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
        solutionLocalView,
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
