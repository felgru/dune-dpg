// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONSTRAINEDLOCALINDEXSET_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONSTRAINEDLOCALINDEXSET_HH

#include <dune/functions/functionspacebases/subspacelocalview.hh>


namespace Dune {
namespace Functions {



template<class LV, class NIS>
class ConstrainedLocalIndexSet
{
public:
  using LocalView = LV;
  using NodeIndexSet = NIS;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = typename NodeIndexSet::MultiIndex;
  using ConstraintWeights = typename NodeIndexSet::ConstraintWeights;
  using size_type = std::size_t;


  ConstrainedLocalIndexSet(const NodeIndexSet& nodeIndexSet) :
    nodeIndexSet_(nodeIndexSet)
  {}

  ConstrainedLocalIndexSet(NodeIndexSet&& nodeIndexSet) :
    nodeIndexSet_(nodeIndexSet)
  {}

  /** \brief Bind the index set to a LocalView
   */
  void bind(const LocalView& localView)
  {
    localView_ = &localView;
    nodeIndexSet_.bind(localView_->tree());
  }

  /** \brief Bind the index set to a SubspaceLocalView
   */
  template<class TreePath>
  void bind(const SubspaceLocalView<LocalView, TreePath>& subspaceLocalView)
  {
    bind(subspaceLocalView.rootLocalView());
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    localView_ = nullptr;
    nodeIndexSet_.unbind();
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    return nodeIndexSet_.size();
  }

  const std::vector<MultiIndex>& indicesLocalGlobal() const
  {
    return nodeIndexSet_.indicesLocalGlobal();
  }

  size_type constraintsSize() const
  {
    return nodeIndexSet_.constraintsSize();
  }

  size_type constraintOffset(size_type i) const
  {
    return nodeIndexSet_.constraintOffset(i);
  }

  const ConstraintWeights& constraintWeights(size_type i) const
  {
    return nodeIndexSet_.constraintWeights(i);
  }

  /** \brief Return the local view that we are attached to
   */
  const LocalView& localView() const
  {
    return *localView_;
  }

protected:

  const LocalView* localView_;

  NodeIndexSet nodeIndexSet_;
};



} // end namespace Functions
} // end namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONSTRAINEDLOCALINDEXSET_HH
