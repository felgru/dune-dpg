#ifndef DUNE_DPG_BOUNDARY_TOOLS
#define DUNE_DPG_BOUNDARY_TOOLS

#include <iostream>

#include <vector>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dune/istl/io.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discretescalarglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>



namespace Dune {

  class BoundaryTools
  {
  public:
    template <class FEBasis>
    void boundaryTreatmentInflow(
                  const FEBasis& ,
                  std::vector<bool>& ,
                  const FieldVector<double,2>&
                  );
  private:
    std::vector<int> getVertexOfIntersection(
                        int
                        );


  };


  //*******************************************************************
  template <class FEBasis >
  void BoundaryTools::boundaryTreatmentInflow(
                        const FEBasis& feBasis,
                        std::vector<bool>& dirichletNodes,
                        const FieldVector<double,2>& beta
                        )
  {
    const int dim = FEBasis::GridView::dimension;

    // Get the grid view from the finite element basis
    typedef typename FEBasis::GridView GridView;
    GridView gridView = feBasis.gridView();

    const int dofs = feBasis.indexSet().size();
    //std::cout<< "dofs=" << dofs << std::endl<< std::endl;

    dirichletNodes.resize(dofs);
    std::vector<int> dirichletNodesInt(dofs,0);

    //We get the local view
    auto localView = feBasis.localView();
    auto localIndexSet = feBasis.indexSet().localIndexSet();


    // A loop over all elements of the grid
    for(const auto& e : elements(gridView))
    {
      // std::cout << std::endl << "NEW ELEMENT" << std::endl<< std::endl;

      // We bind the local view to the current element
      localView.bind(e);
      const auto& localFEM = localView.tree().finiteElement();

      // We bind the local index set
      localIndexSet.bind(localView);

      // dofs in the current finite element
      int dofsLocal = localFEM.localCoefficients().size();

      int nFace = 3;  // todo: grab this generically!
      int nVertex = nFace;// todo: grab this generically!

      // For every vertex, we have to see whether it is on the inflow boundary
      // If vertex i is on the inflow boundary, we will have vertexOnInflowBoundary[i] >0
      std::vector<int> vertexOnInflowBoundary(nVertex,0);

      // for all intersections, we see which one lies on the inflow boundary
      // if intersection i lies on the inflow boundary, then faceOnInflowBoundary[i]=true
      // we will assume that an intersection is simply a face for us


      std::vector<int> faceOnInflowBoundary(nFace,0);

      for (auto&& intersection : intersections(gridView, e))
      {
        // std::cout << "****** New intersection" << std::endl<< std::endl;

        //Local index of the intersection
        int indexIntersection = intersection.indexInInside();

        //outer normal vector in the center of the face
        const FieldVector<double,dim>& centerOuterNormal =
               intersection.centerUnitOuterNormal();

        // n.beta
        double scalarProd = centerOuterNormal * beta;

        // We see whether we are on the inflow boundary
        double tolerance = -1e-8*beta.two_norm();
        bool isOnInflowBoundary = (scalarProd<tolerance)
                                  && intersection.boundary();

        // We store this information in faceOnInflowBoundary
        faceOnInflowBoundary[indexIntersection] = isOnInflowBoundary;

        // std::cout << "indexIntersection=" << indexIntersection << std::endl;
        // std::cout << "isOnInflowBoundary=" << isOnInflowBoundary << std::endl<< std::endl;

        //if the intersection is on the inflow boundary, we have to update what are the local
        // verteces that are also on the inflow boundary
        if(isOnInflowBoundary)
        {
          // We see what are the vertices associated to the current intersection (assumed to be a face)
          // TODO: This is hard coded and works only for triangles so it would be good to implement it better
          std::vector<int> vertexOfIntersection = this->getVertexOfIntersection(indexIntersection);

          vertexOnInflowBoundary[ vertexOfIntersection[0] ] += 1;
          vertexOnInflowBoundary[ vertexOfIntersection[1] ] += 1;

          // std::cout << "vertexOfIntersection[" << 0 << "]=" << vertexOfIntersection[0] << std::endl;
          // std::cout << "vertexOfIntersection[" << 1 << "]=" << vertexOfIntersection[1] << std::endl;
        }

        // std::cout << "scalarProd=" << scalarProd << "; isOnInflowBoundary=" << isOnInflowBoundary << std::endl;
        // std::cout << "centerOuterNormal[" << 0 << "]=" << centerOuterNormal[0] << std::endl;
        // std::cout << "centerOuterNormal[" << 1 << "]=" << centerOuterNormal[1] << std::endl << std::endl;

      }



      //  for(int i=0;i<vertexOnInflowBoundary.size();i++){
      // std::cout << "vertexOnInflowBoundary[" << i << "]=" << vertexOnInflowBoundary[i] << std::endl;

      //  }
      //  for(int i=0;i<faceOnInflowBoundary.size();i++){
      // std::cout << "faceOnInflowBoundary[" << i << "]=" << faceOnInflowBoundary[i] << std::endl;

      //  }

      //For each dof, we check whether it belongs to the inflow boundary
      for(int i=0; i<dofsLocal; i++)
      {

        int dofOnInflowBoundary = 0;

        //localkey of dof i
        const auto& dofLocalKey = localFEM.localCoefficients().localKey(i);

        // Codimension and subentity index of the current dof
        int dofCodim = dofLocalKey.codim();
        int dofIndex = dofLocalKey.subEntity();

        if(dofCodim == 1)//the dof belongs to a face
        {
           dofOnInflowBoundary += faceOnInflowBoundary[dofIndex];
        }
        if(dofCodim == 2)//the dof belongs to a vertex
        {
          dofOnInflowBoundary += vertexOnInflowBoundary[dofIndex];
        }

        // std::cout <<"dofs: local index= " << i
        //       << " ; global index = " << localIndexSet.index(i)[0]
        //       << " ; codim=" << dofCodim
        //       << " ; subentityIndex=" << dofIndex
        //       << std::endl;

        dirichletNodesInt[ localIndexSet.index(i)[0] ] += dofOnInflowBoundary;


      } // end dof


    } // end element e

    for(int i=0; i<dirichletNodes.size(); i++)
    {
      dirichletNodes[i] = (dirichletNodesInt[i]>0);
    }

    return;
  }

//*************************************
  std::vector<int> BoundaryTools::getVertexOfIntersection(
                              int indexIntersection
                              )
  {
    std::vector<int> indexVertex(2,0);

    if(indexIntersection==0)
    {
      indexVertex[0]=0;
      indexVertex[1]=1;
    }
    if(indexIntersection==1)
    {
      indexVertex[0]=0;
      indexVertex[1]=2;
    }
    if(indexIntersection==2)
    {
      indexVertex[0]=1;
      indexVertex[1]=2;
    }

    return indexVertex;

  }



} // end namespace Dune

#endif // DUNE_DPG_BOUNDARY_TOOLS
