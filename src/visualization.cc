#include <config.h>

#include <string>
#include <vector>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/istl/bvector.hh>

#include <dune/localfunctions/lagrange/q1.hh>
#include <dune/localfunctions/lagrange/qk.hh>

#include <dune/localfunctions/test/test-localfe.hh>

#include <dune/dpg/visualization.hh>

using namespace Dune;


int main()
{

    // ////////////////////////////////
    //   Generate the grid
    // ////////////////////////////////

    constexpr int dim = 2;
    typedef YaspGrid<dim> GridType;
    const FieldVector<double,dim> l(1);
    const GridType grid(l, {{1,1}});
    const GridType fineGrid(l, {{10,10}});

    // ///////////////////////////////////////////////////////
    //   Stiffness matrix and right hand side vector
    // ///////////////////////////////////////////////////////

    typedef BlockVector<FieldVector<double,1> > VectorType;

    // ///////////////////////////////////////////////////////
    //   Evaluate the basis functions
    // ///////////////////////////////////////////////////////

    QkLocalFiniteElement<double,double,dim,2> localFiniteElement;
    std::vector<VectorType> basisFunctions;
    evaluateBasisFunctions(localFiniteElement,
                           grid.leafGridView(),
                           fineGrid.leafGridView(),
                           basisFunctions);

    testFE(localFiniteElement);

    for(unsigned int i = 0; i < basisFunctions.size(); i++) {
      // Output result
      const std::string description = "basis function " + std::to_string(i);
      const std::string filename    = "basis_function_" + std::to_string(i);
      VTKWriter<GridType::LeafGridView> vtkWriter(fineGrid.leafGridView());
      vtkWriter.addVertexData(basisFunctions[i], description.c_str());
      vtkWriter.write(filename.c_str());
    }
}
