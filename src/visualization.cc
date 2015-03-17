#include <config.h>

#include <vector>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/istl/bcrsmatrix.hh>

#include <dune/localfunctions/lagrange/q1.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/localfunctions/lagrange/qkdiscontinuous.hh>

#include <dune/localfunctions/test/test-localfe.hh>

#include <dune/dpg/visualization.hh>

using namespace Dune;


int main (int argc, char *argv[]) try
{

    // ////////////////////////////////
    //   Generate the grid
    // ////////////////////////////////

    const int dim = 2;
    typedef YaspGrid<dim> GridType;
    FieldVector<double,dim> l(1);
    GridType grid(l, {{1,1}});
    GridType fineGrid(l, {{10,10}});

    // ///////////////////////////////////////////////////////
    //   Stiffness matrix and right hand side vector
    // ///////////////////////////////////////////////////////

    typedef BlockVector<FieldVector<double,1> > VectorType;

    // ///////////////////////////////////////////////////////
    //   Evaluate the basis functions
    // ///////////////////////////////////////////////////////

    QkDiscontinuousLocalFiniteElement<double,double,dim,2> localFiniteElement;
    std::vector<VectorType> basisFunctions;
    evaluateBasisFunctions(localFiniteElement,
                           grid.leafGridView(),
                           fineGrid.leafGridView(),
                           basisFunctions);

    testFE(localFiniteElement);

    // Output result
    VTKWriter<GridType::LeafGridView> vtkWriter0(fineGrid.leafGridView());
    vtkWriter0.addVertexData(basisFunctions[0], "basis function 0");
    vtkWriter0.write("basis_function_0");
    // Output result
    VTKWriter<GridType::LeafGridView> vtkWriter1(fineGrid.leafGridView());
    vtkWriter1.addVertexData(basisFunctions[1], "basis function 1");
    vtkWriter1.write("basis_function_1");
    // Output result
    VTKWriter<GridType::LeafGridView> vtkWriter2(fineGrid.leafGridView());
    vtkWriter2.addVertexData(basisFunctions[2], "basis function 2");
    vtkWriter2.write("basis_function_2");
    // Output result
    VTKWriter<GridType::LeafGridView> vtkWriter3(fineGrid.leafGridView());
    vtkWriter3.addVertexData(basisFunctions[3], "basis function 3");
    vtkWriter3.write("basis_function_3");
    // Output result
    VTKWriter<GridType::LeafGridView> vtkWriter4(fineGrid.leafGridView());
    vtkWriter4.addVertexData(basisFunctions[4], "basis function 4");
    vtkWriter4.write("basis_function_4");
    // Output result
    VTKWriter<GridType::LeafGridView> vtkWriter5(fineGrid.leafGridView());
    vtkWriter5.addVertexData(basisFunctions[5], "basis function 5");
    vtkWriter5.write("basis_function_5");
    // Output result
    VTKWriter<GridType::LeafGridView> vtkWriter6(fineGrid.leafGridView());
    vtkWriter6.addVertexData(basisFunctions[6], "basis function 6");
    vtkWriter6.write("basis_function_6");
    // Output result
    VTKWriter<GridType::LeafGridView> vtkWriter7(fineGrid.leafGridView());
    vtkWriter7.addVertexData(basisFunctions[7], "basis function 7");
    vtkWriter7.write("basis_function_7");
    // Output result
    VTKWriter<GridType::LeafGridView> vtkWriter8(fineGrid.leafGridView());
    vtkWriter8.addVertexData(basisFunctions[8], "basis function 8");
    vtkWriter8.write("basis_function_8");
}
// Error handling
catch (Exception e) {
    std::cout << e << std::endl;
}
