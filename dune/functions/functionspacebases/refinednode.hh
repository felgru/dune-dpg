// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDNODE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDNODE_HH

#include <dune/common/hybridutilities.hh>
#include <dune/common/math.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/functions/functionspacebases/referencerefinementfactory.hh>
#include <dune/typetree/traversal.hh>
#include <dune/typetree/visitor.hh>

namespace Dune {
namespace Functions {

template<typename Element, typename ctype, int dim, int level>
class RefinedNode
{
  using RefinementCache = ReferenceRefinementCache<ctype, dim, level>;
  using GridType = typename RefinementCache::GridType;

public:

  using RefinementGrid = GridType;
  using RefinementGridView = typename RefinementGrid::LeafGridView;
  using SubElement = typename RefinementGrid::template Codim<0>::Entity;

  RefinedNode() :
    element_(nullptr),
    subElement_(nullptr)
  {}

  RefinementGridView refinedReferenceElementGridView() const
  {
    return refinementCache_.get(element_->type()).leafGridView();
  }

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *element_;
  }

  const SubElement& subElement() const
  {
    return *subElement_;
  }

protected:

  static RefinementCache refinementCache_;
  const Element* element_;
  const SubElement* subElement_;
};

template<typename Element, typename ctype, int dim, int level>
ReferenceRefinementCache<ctype, dim, level>
    RefinedNode<Element, ctype, dim, level>::refinementCache_{};

template<int dim, int level, int k>
struct DGRefinedPreBasisConstants {};

template<int level, int k>
struct DGRefinedPreBasisConstants<2, level, k>
{
  const static int dofsPerSubEdge        = k+1;
  const static int dofsPerSubTriangle    = (k+1)*(k+2)/2;
  const static int dofsPerSubQuad        = (k+1)*(k+1);

  const static int numberOfSubEdges      = power(2, level);
  const static int numberOfSubTriangles  = power(4, level);
  const static int numberOfSubQuads      = power(4, level);
};

class ResetSubElementsVisitor
  : public TypeTree::TreeVisitor
  , public TypeTree::DynamicTraversal
{
  template<typename Node_>
  using hasResetSubElements
      = decltype(std::declval<Node_>().resetSubElements());

public:
  template<typename Node, typename TreePath>
  void pre(Node& node, TreePath treePath)
  {}

  template<typename Node, typename TreePath>
  void post(Node& node, TreePath treePath)
  {}

  template<typename Node, typename TreePath>
  void leaf(Node& node, TreePath treePath)
  {
    if constexpr (Std::is_detected<hasResetSubElements, Node>{})
      node.resetSubElements();
  }
};

template<typename Tree>
void resetSubElementsOfTree(Tree& tree)
{
  ResetSubElementsVisitor visitor{};
  TypeTree::applyToTree(tree,visitor);
}

template<typename Entity>
struct BindToSubElementVisitor
  : public TypeTree::TreeVisitor
  , public TypeTree::DynamicTraversal
{

  template<typename Node, typename TreePath>
  void pre(Node& node, TreePath treePath)
  {}

  template<typename Node, typename TreePath>
  void post(Node& node, TreePath treePath)
  {}

  template<typename Node, typename TreePath>
  void leaf(Node& node, TreePath treePath)
  {
    node.bindSubElement(subEntity_);
  }

  BindToSubElementVisitor(const Entity& subEntity)
    : subEntity_(subEntity)
  {}

  const Entity& subEntity_;

};

template<typename Tree, typename Entity>
void bindTreeToSubElement(Tree& tree, const Entity& subEntity)
{
  BindToSubElementVisitor<Entity> visitor(subEntity);
  TypeTree::applyToTree(tree,visitor);
}

}} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDNODE_HH
