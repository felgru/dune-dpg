#ifndef DUNE_DPG_VISULAIZATION_HH
#define DUNE_DPG_VISULAIZATION_HH

#include <vector>

#include <dune/istl/bvector.hh>


namespace Dune
{

/** \brief Evaluate all basis functions on the given coarse grid view on the fine grid view*/
template <class LocalFiniteElement, class GridView>
void evaluateBasisFunctions(
          const LocalFiniteElement& localFiniteElement,
          const GridView& coarseGridView,
          const GridView& fineGridView,
          std::vector<BlockVector<FieldVector<double,1> > >& basisEvaluations)
{
    constexpr int dim = GridView::dimension;

    // The index set gives you indices for each element, edge, face, vertex, etc.
    const typename GridView::IndexSet& indexSet = fineGridView.indexSet();

    // This only works for cube grids
    assert(localFiniteElement.type()
           == coarseGridView.template begin<0>()->type());

    // set basisEvaluations to correct length
    basisEvaluations.resize(localFiniteElement.size());
    for(auto& v : basisEvaluations)
    {
        v.resize(fineGridView.size(dim));
        // Set all entries to zero
        v = 0;
    }

    std::vector<FieldVector<double,1>> evaluation;
    evaluation.reserve(localFiniteElement.size());

    for(const auto& v : vertices(fineGridView)) {
        // The shape functions on the reference element evaluated on the
        // finer grid
        const FieldVector<double,dim> vPos = v.geometry().corner(0);
        localFiniteElement.localBasis().evaluateFunction(vPos, evaluation);

        for(size_t i=0; i<evaluation.size(); i++) {
            // The global index of the vertex v
            const auto vertexIndex = indexSet.index(v);

            for(size_t j=0; j<evaluation.size(); j++)
                basisEvaluations[i][vertexIndex] = evaluation[i];
        }
    }
}

} // end namespace Dune

#endif // DUNE_DPG_VISULAIZATION_HH
