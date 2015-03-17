#include <vector>

#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/function.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/common/intersection.hh> //TODO necessary?

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dune/common/std/final.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

namespace Dune {

// Compute the stiffness matrix for a single element
template <class LocalViewTrace, class LocalViewInterior, class LocalViewTest, class GridView, class MatrixType>
void getLocalMatrix(const LocalViewTrace& localViewTrace,
                    const LocalViewInterior& localViewInterior,
                    const LocalViewTest& localViewTest,
                    const GridView& gridView,
                    MatrixType& elementMatrix)
{
  FieldVector<double, 2> beta={1,1};
  double c=0;

  // Get the grid element from the local FE basis view
  typedef typename LocalViewTest::Element Element;
  const Element& element = localViewTest.element();

  const int dim = Element::dimension;
  auto geometry = element.geometry();

  // Get set of shape functions for this element
  const auto& localFiniteElementTest = localViewTest.tree().finiteElement();
  const auto& localFiniteElementTrace = localViewTrace.tree().finiteElement();
  const auto& localFiniteElementInterior = localViewInterior.tree().finiteElement();

  const int nTest(localFiniteElementTest.localBasis().size());
  const int nTrace(localFiniteElementTrace.localBasis().size());
  const int nInterior(localFiniteElementInterior.localBasis().size());

  // Set all matrix entries to zero
  elementMatrix.setSize(nTest, nTest+nTrace+nInterior);
  elementMatrix = 0;      // fills the entire matrix with zeroes

  ////////////////////////////
  // Assemble interior terms
  ////////////////////////////

  // Get a quadrature rule
  int order = 2*(dim*localFiniteElementTest.localBasis().order()-1);  //TODO!!!!!!!!!
  const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (size_t pt=0; pt < quad.size(); pt++) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The transposed inverse Jacobian of the map from the reference element to the element
    const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPos);

    ////////////////////////////
    // Test Functions:
    ////////////////////////////
    // values of the shape functions
    std::vector<FieldVector<double,1> > valuesTest;
    localFiniteElementTest.localBasis().evaluateFunction(quadPos, valuesTest);

    // The gradients of the shape functions on the reference element
    std::vector<FieldMatrix<double,1,dim> > referenceGradientsTest;
    localFiniteElementTest.localBasis().evaluateJacobian(quadPos, referenceGradientsTest);

    // Compute the shape function gradients on the real element
    std::vector<FieldVector<double,dim> > gradientsTest(referenceGradientsTest.size());
    for (size_t i=0; i<gradientsTest.size(); i++)
      jacobian.mv(referenceGradientsTest[i][0], gradientsTest[i]);

    ////////////////////////////
    // Interior Trial Functions
    ////////////////////////////
    // values of the shape functions
    std::vector<FieldVector<double,1> > valuesInterior;
    localFiniteElementInterior.localBasis().evaluateFunction(quadPos, valuesInterior);

    // Compute the actual matrix entries
    for (size_t i=0; i<nTest; i++)
    {
      for (size_t j=0; j<nTest; j++ )
      {
        elementMatrix[i][j] += (valuesTest[i] * valuesTest[j]) * quad[pt].weight() * integrationElement;
        elementMatrix[i][j] += (beta*gradientsTest[i]) * (beta*gradientsTest[j]) * quad[pt].weight() * integrationElement;
      }
      for (size_t j=0; j<nInterior; j++)
      {
        elementMatrix[i][j+nTest+nTrace] += (valuesTest[i] * valuesInterior[j])* c * quad[pt].weight() * integrationElement;
        elementMatrix[i][j+nTest+nTrace] += (beta*gradientsTest[i])*(-1.0) * valuesInterior[j] * quad[pt].weight() * integrationElement;
      }
    }
  }


  ////////////////////////////
  // Assemble boundary terms
  ////////////////////////////

  for (auto&& intersection : intersections(gridView, element))
  {
    const QuadratureRule<double, dim-1>& quadFace = QuadratureRules<double, dim-1>::rule(intersection.type(), order);
    // Loop over all quadrature points
    for (size_t pt=0; pt < quadFace.size(); pt++) {

    // Position of the current quadrature point in the reference element (face!)
    const FieldVector<double,dim-1>& quadFacePos = quadFace[pt].position();

    // The multiplicative factor in the integral transformation formula moltiplied with outer normal
    const FieldVector<double,dim>& integrationOuterNormal = intersection.integrationOuterNormal(quadFacePos);

                // position of the quadrature point within the element
    const FieldVector<double,dim> elementQuadPos = intersection.geometryInInside().global(quadFacePos);


    ////////////////////////////
    // Test Functions:
    ////////////////////////////
    // values of the shape functions
    std::vector<FieldVector<double,1> > valuesTest;
    localFiniteElementTest.localBasis().evaluateFunction(elementQuadPos, valuesTest);

    ////////////////////////////
    // Trace Functions
    ////////////////////////////
    // values of the shape functions
    std::vector<FieldVector<double,1> > valuesTrace;
    localFiniteElementTrace.localBasis().evaluateFunction(elementQuadPos, valuesTrace);

    // Compute the actual matrix entries
    for (size_t i=0; i<nTest; i++)
    {
      for (size_t j=0; j<nTrace; j++ )
      {
        elementMatrix[i][j+nTest] += ((beta*integrationOuterNormal) * valuesTest[i] * valuesTrace[j]) * quadFace[pt].weight();
      }
    }
    }
  }

}



// Compute the source term for a single element
template <class LocalViewTest, class LocalVolumeTerm>
void getVolumeTerm( const LocalViewTest& localViewTest,
                    BlockVector<FieldVector<double,1> >& localRhs,
                    //const Dune::VirtualFunction<FieldVector<double,LocalViewTest::Element::dimension>, double>* volumeTerm)
                    LocalVolumeTerm&& localVolumeTerm)
{
  // Get the grid element from the local FE basis view
  typedef typename LocalViewTest::Element Element;
  const Element& element = localViewTest.element();

  const int dim = Element::dimension;

  // Get set of shape functions for this element
  const auto& localFiniteElementTest = localViewTest.tree().finiteElement();

  // Set all entries to zero
  localRhs.resize(localFiniteElementTest.localBasis().size());
  localRhs = 0;

  // A quadrature rule
  int order = dim*localFiniteElementTest.localBasis().order(); //TODO!!!!!!
  const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element.type(), order);


  // Loop over all quadrature points
  for ( size_t pt=0; pt < quad.size(); pt++ ) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = element.geometry().integrationElement(quadPos);

    double functionValue = localVolumeTerm(quadPos);

    // Evaluate all shape function values at this point
    std::vector<FieldVector<double,1> > shapeFunctionValues;
    localFiniteElementTest.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

    // Actually compute the vector entries
    for (size_t i=0; i<localRhs.size(); i++)
      localRhs[i] += shapeFunctionValues[i] * functionValue * quad[pt].weight() * integrationElement;

  }

}

// Get the occupation pattern of the stiffness matrix TODO!!!!!!!!
template <class FEBasisTrace, class FEBasisInterior, class FEBasisTest>
void getOccupationPattern(const FEBasisTrace& feBasisTrace,
                            const FEBasisInterior& feBasisInterior,
                            const FEBasisTest& feBasisTest,
                            MatrixIndexSet& nb)
{
  // Total number of degrees of freedom
  auto basisIndexSetTrace = feBasisTrace.indexSet();
  auto basisIndexSetInterior = feBasisInterior.indexSet();
  auto basisIndexSetTest = feBasisTest.indexSet();
  auto n = basisIndexSetTrace.size()+basisIndexSetInterior.size()+basisIndexSetTest.size();

  nb.resize(n, n);

  // A view on the FE basis on a single element
  typename FEBasisTest::LocalView localViewTest(&feBasisTest);
  typename FEBasisTrace::LocalView localViewTrace(&feBasisTrace);
  typename FEBasisInterior::LocalView localViewInterior(&feBasisInterior);

  auto localIndexSetTrace = basisIndexSetTrace.localIndexSet();
  auto localIndexSetInterior = basisIndexSetInterior.localIndexSet();
  auto localIndexSetTest = basisIndexSetTest.localIndexSet();

  // Loop over all leaf elements
  auto it    = feBasisTest.gridView().template begin<0>();
  auto endIt = feBasisTest.gridView().template end<0>  ();
  auto itTrace    = feBasisTrace.gridView().template begin<0>();
  auto itInterior    = feBasisInterior.gridView().template begin<0>();

  int offsetTrace(feBasisTest.indexSet().size());
  int offsetInterior(offsetTrace+feBasisTrace.indexSet().size());

  for (; it!=endIt; ++it, ++itTrace, ++itInterior)
  {
    // Bind the local FE basis view to the current element
    localViewTest.bind(*it);
    localViewTrace.bind(*itTrace);
    localViewInterior.bind(*itInterior);

    localIndexSetTrace.bind(localViewTrace);
    localIndexSetInterior.bind(localViewInterior);
    localIndexSetTest.bind(localViewTest);

        // There is a matrix entry a_ij if the i-th and j-th vertex are connected in the grid
        for (size_t i=0; i<localViewTest.tree().size(); i++) {

            for (size_t j=0; j<localViewTest.tree().size(); j++) {

                auto iIdx = localIndexSetTest.index(i)[0];
                auto jIdx = localIndexSetTest.index(j)[0];

                // Add a nonzero entry to the matrix
                nb.add(iIdx, jIdx);

            }
            for (size_t j=0; j<localViewTrace.tree().size(); j++) {

                auto iIdx = localIndexSetTest.index(i)[0];
                auto jIdx = localIndexSetTrace.index(j)[0];

                // Add a nonzero entry to the matrix
                nb.add(iIdx, jIdx+offsetTrace);
                nb.add(jIdx+offsetTrace, iIdx);

            }
            for (size_t j=0; j<localViewInterior.tree().size(); j++) {

                auto iIdx = localIndexSetTest.index(i)[0];
                auto jIdx = localIndexSetInterior.index(j)[0];

                // Add a nonzero entry to the matrix
                nb.add(iIdx, jIdx+offsetInterior);
                nb.add(jIdx+offsetInterior, iIdx);
            }
        }
    }
    for (size_t i=0; i<n; i++)
    {
      nb.add(i, i);
    }

}


/** \brief Assemble the Laplace stiffness matrix on the given grid view */
template <class FEBasisTrace, class FEBasisInterior, class FEBasisTest, class VolumeTerm>
void assembleSystem(const FEBasisTrace& feBasisTrace,
                    const FEBasisInterior& feBasisInterior,
                    const FEBasisTest& feBasisTest,
                    BCRSMatrix<FieldMatrix<double,1,1> >& matrix,
                    BlockVector<FieldVector<double,1> >& rhs,
                    VolumeTerm&& volumeTerm)
{
  // Get the grid view from the finite element basis
  typedef typename FEBasisTest::GridView GridView;
  GridView gridView = feBasisTest.gridView();

  auto localVolumeTerm = localFunction(Functions::makeGridViewFunction(volumeTerm, gridView));

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  MatrixIndexSet occupationPattern;
  getOccupationPattern<FEBasisTrace, FEBasisInterior, FEBasisTest>(feBasisTrace, feBasisInterior, feBasisTest, occupationPattern);

  // ... and give it the occupation pattern we want.
  occupationPattern.exportIdx(matrix);

  // set rhs to correct length -- the total number of basis vectors in the basis
  auto basisIndexSetTrace = feBasisTrace.indexSet();
  auto basisIndexSetInterior = feBasisInterior.indexSet();
  auto basisIndexSetTest = feBasisTest.indexSet();
  rhs.resize(basisIndexSetTrace.size()+basisIndexSetInterior.size()+basisIndexSetTest.size());

  // Set all entries to zero
  matrix = 0;
  rhs = 0;

  // A view on the FE basis on a single element
  auto localViewTrace = feBasisTrace.localView();
  auto localViewInterior = feBasisInterior.localView();
  auto localViewTest = feBasisTest.localView();

  auto localIndexSetTrace = basisIndexSetTrace.localIndexSet();
  auto localIndexSetInterior = basisIndexSetInterior.localIndexSet();
  auto localIndexSetTest = basisIndexSetTest.localIndexSet();

  // A loop over all elements of the grid
  for(const auto& e : elements(gridView)) {

    // Bind the local FE basis view to the current element
    localViewTrace.bind(e);
    localViewInterior.bind(e);
    localViewTest.bind(e);
    localIndexSetTrace.bind(localViewTrace);
    localIndexSetInterior.bind(localViewInterior);
    localIndexSetTest.bind(localViewTest);

    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    getLocalMatrix(localViewTrace, localViewInterior, localViewTest, gridView, elementMatrix);

    const int nTest(localViewTest.tree().finiteElement().localBasis().size());
    const int nTrace(localViewTrace.tree().finiteElement().localBasis().size());
    const int nInterior(localViewInterior.tree().finiteElement().localBasis().size());

    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i=0; i<nTest; i++)
    {
      // The global index of the i-th local degree of freedom of the element 'it'
      auto row = localIndexSetTest.index(i)[0];

      for (size_t j=0; j<nTest; j++ )
      {
        // The global index of the j-th local degree of freedom of the element 'it'
        auto col = localIndexSetTest.index(j)[0];
        matrix[row][col] += elementMatrix[i][j];
      }
      for (size_t j=0; j<nTrace; j++ )
      {
        // The global index of the j-th local degree of freedom of the element 'it'
        auto col = localIndexSetTrace.index(j)[0]
                    +feBasisTest.indexSet().size();
        matrix[row][col] += elementMatrix[i][j+nTest];
        matrix[col][row] += elementMatrix[i][j+nTest];
      }
      for (size_t j=0; j<nInterior; j++ )
      {
        // The global index of the j-th local degree of freedom of the element 'it'
        auto col = localIndexSetInterior.index(j)[0]
                    +feBasisTest.indexSet().size()
                    +feBasisTrace.indexSet().size();
        matrix[row][col] += elementMatrix[i][j+nTest+nTrace];
        matrix[col][row] += elementMatrix[i][j+nTest+nTrace];
      }
    }

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> > localRhs;
    localVolumeTerm.bind(e);
    getVolumeTerm(localViewTest, localRhs, localVolumeTerm);

     for (size_t i=0; i<localRhs.size(); i++) {

      // The global index of the i-th vertex of the element 'it'
      auto row = localIndexSetTest.index(i)[0];
      rhs[row] += localRhs[i];

    }

  }

}

} // end namespace Dune
