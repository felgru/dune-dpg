// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDNODE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDNODE_HH

#include <dune/common/power.hh>
#include <dune/functions/functionspacebases/referencerefinementfactory.hh>

namespace Dune {
namespace Functions {

template<typename Element, typename ctype, int dim, int level>
class RefinedNode
{
  using RefinementCache = ReferenceRefinementCache<ctype, dim, level>;
  using GridType = typename RefinementCache::GridType;

public:

  using RefinementGrid = GridType;

  RefinedNode() :
    element_(nullptr)
  {}

  const RefinementGrid& refinedReferenceElement() const
  {
    return refinementCache_.get(element_->type());
  }

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *element_;
  }

protected:

  static RefinementCache refinementCache_;
  const Element* element_;
};

template<typename Element, typename ctype, int dim, int level>
ReferenceRefinementCache<ctype, dim, level>
    RefinedNode<Element, ctype, dim, level>::refinementCache_{};

template<int dim, int level, int k>
struct DGRefinedNodeFactoryConstants {};

template<int level, int k>
struct DGRefinedNodeFactoryConstants<2, level, k>
{
  const static int dofsPerSubEdge        = k+1;
  const static int dofsPerSubTriangle    = (k+1)*(k+2)/2;
  const static int dofsPerSubQuad        = (k+1)*(k+1);

  const static int numberOfSubEdges        = StaticPower<2,level>::power;
  const static int numberOfSubTriangles    = StaticPower<4,level>::power;
  const static int numberOfSubQuads        = StaticPower<4,level>::power;
};

}} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REFINEDNODE_HH
