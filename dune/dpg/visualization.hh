#ifndef DUNE_DPG_VISULAIZATION_HH
#define DUNE_DPG_VISULAIZATION_HH

#include <vector>

#include <dune/istl/bcrsmatrix.hh>


namespace Dune
{

/** \brief Evaluate all basis functions on the given coarse grid view on the fine grid view*/
template <class LocalFiniteElement, class GridView>
void evaluateBasisFunctions(
          const LocalFiniteElement& localFiniteElement,
          const GridView& gridView,
          const GridView& fineGridView,
          std::vector<BlockVector<FieldVector<double,1> > >& basisEvaluations)
{
    const int dim = GridView::dimension;

    // The index set gives you indices for each element, edge, face, vertex, etc.
    const typename GridView::IndexSet& indexSet = fineGridView.indexSet();

    // The element of the coarse 1-cell grid
    auto elemIt = gridView.template begin<0>();

    auto geometry = elemIt->geometry();

    assert(localFiniteElement.type() == elemIt->type());  // This only works for cube grids

    // set basisEvaluations to correct length
    basisEvaluations.resize(localFiniteElement.size());
    for(auto it    = basisEvaluations.begin(),
             endIt = basisEvaluations.end()  ;
        it != endIt; ++it)
    {
        it->resize(fineGridView.size(dim));
        // Set all entries to zero
        *it = 0;
    }

    // A loop over all vertices of the fine grid
    auto it    = fineGridView.template begin<dim>();
    auto endIt = fineGridView.template end<dim>  ();

    for( ; it != endIt; ++it ) {

        // The shape functions on the reference element evaluated on the
        // finer grid
        std::vector<FieldVector<double,1> > evaluation;
        FieldVector<double,dim> vPos = it->geometry().corner(0);
        localFiniteElement.localBasis().evaluateFunction(vPos, evaluation);

        //
        for(size_t i=0; i<evaluation.size(); i++) {

            // The global index of the vertex given by 'it'
            auto vertexIndex = indexSet.index(*it);

            for(size_t j=0; j<evaluation.size(); j++)
                basisEvaluations[i][vertexIndex] = evaluation[i];
        }

    }

}

} // end namespace Dune

#endif // DUNE_DPG_VISULAIZATION_HH
