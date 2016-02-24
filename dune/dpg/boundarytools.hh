#ifndef DUNE_DPG_BOUNDARY_TOOLS
#define DUNE_DPG_BOUNDARY_TOOLS

#include <iostream>

#include <vector>

#include <dune/common/fvector.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dune/istl/io.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>



namespace Dune {

  template<class Function>
  class BoundaryCondition{

  Function g_;

  public:

    using DomainType = FieldVector<double, 2>;
    using RangeType  = FieldVector<double, 1>;

    BoundaryCondition(Function g) : g_(g) {};

    void evaluate(
      const DomainType& x,
      RangeType& y) const; /// Remark: this signature assumes that we have a 2D scalar problem

  };

  class BoundaryTools
  {

  public:
    template <class FEBasis,class Direction>
    static void getInflowBoundaryMask(
                  const FEBasis& ,
                  std::vector<bool>& ,
                  const Direction&
                  );
    template <class FEBasis, class Function>
    static void getInflowBoundaryValue(
                  const FEBasis& ,
                  std::vector<double>& ,
                  Function&
                  );

  private:
    static std::vector<unsigned int> getVertexOfIntersection(
                        unsigned int ,
                        GeometryType
                        );

  };

  //*******************************************************************
  template<class Function>
  void BoundaryCondition<Function>::evaluate(
                                            const DomainType& x,
                                            RangeType& y
                                            ) const
  {
    y = std::get<0>(g_)(x);
  }

  //*******************************************************************
  /**
 * \brief Writes in the vector dirichletNodes whether a degree of freedom is in the boundary (value 1) or not (value 0).
 *        The degrees of freedom are relative to the finite element basis feBasis.
 *        In the transport problem, the result depends on the direction of propagation beta.
 * \param feBasis        a finite element basis
 * \param dirichletNodes the vector where we store the output
 * \param beta           the direction of propagation
 */
  template <class FEBasis,class Direction>
  void BoundaryTools::getInflowBoundaryMask(
                        const FEBasis& feBasis,
                        std::vector<bool>& dirichletNodes,
                        const Direction& beta
                        )
  {
    const unsigned int dim = FEBasis::GridView::dimension;

    typedef typename FEBasis::GridView GridView;
    GridView gridView = feBasis.gridView();

    const unsigned int dofs = feBasis.size();

    dirichletNodes.resize(dofs);
    std::vector<unsigned int> dirichletNodesInt(dofs,0);

    auto localView = feBasis.localView();
    auto localIndexSet = feBasis.localIndexSet();


    for(const auto& e : elements(gridView))
    {
      localView.bind(e);
      const auto& localFEM = localView.tree().finiteElement();

      localIndexSet.bind(localView);

      // dofs in the current finite element
      const unsigned int dofsLocal = localFEM.localCoefficients().size();

      const unsigned int nFace
          = ReferenceElements<double, dim>::general(e.type()).size(dim-1);
      const unsigned int nVertex
          = ReferenceElements<double, dim>::general(e.type()).size(dim);

      // For every vertex, we have to see whether it is on the inflow boundary.
      // If vertex i is on the inflow boundary, we will have vertexOnInflowBoundary[i] >0.
      std::vector<unsigned int> vertexOnInflowBoundary(nVertex,0);

      // for all intersections, we see which one lies on the inflow boundary
      // if intersection i lies on the inflow boundary, then faceOnInflowBoundary[i]=true
      // we will assume that an intersection is simply a face for us


      std::vector<unsigned int> faceOnInflowBoundary(nFace,0);

      for (auto&& intersection : intersections(gridView, e))
      {
        // Local index of the intersection
        const unsigned int indexIntersection = intersection.indexInInside();

        // outer normal vector in the center of the face
        const FieldVector<double,dim>& centerOuterNormal =
               intersection.centerUnitOuterNormal();

        // n.beta
        const double scalarProd = centerOuterNormal * beta;

        // We see whether we are on the inflow boundary
        const double tolerance = -1e-8*beta.two_norm();
        const bool isOnInflowBoundary = (scalarProd < tolerance)
                                        && intersection.boundary();

        // We store this information in faceOnInflowBoundary
        faceOnInflowBoundary[indexIntersection] = isOnInflowBoundary;

        // if the intersection is on the inflow boundary, we have to update
        // what are the local vertices that are also on the inflow boundary
        if(isOnInflowBoundary)
        {
          // We see what are the vertices associated to the current
          // intersection (assumed to be a face)
          std::vector<unsigned int> vertexOfIntersection
              = getVertexOfIntersection(indexIntersection, e.type());

          vertexOnInflowBoundary[ vertexOfIntersection[0] ] += 1;
          vertexOnInflowBoundary[ vertexOfIntersection[1] ] += 1;
        }
      }

      // For each dof, we check whether it belongs to the inflow boundary
      for(unsigned int i=0; i<dofsLocal; i++)
      {
        unsigned int dofOnInflowBoundary = 0;

        // localkey of dof i
        const auto& dofLocalKey = localFEM.localCoefficients().localKey(i);

        // Codimension and subentity index of the current dof
        const unsigned int dofCodim = dofLocalKey.codim();
        const unsigned int dofIndex = dofLocalKey.subEntity();

        if(dofCodim == 1) //the dof belongs to a face
        {
           dofOnInflowBoundary += faceOnInflowBoundary[dofIndex];
        }
        if(dofCodim == 2) //the dof belongs to a vertex
        {
          dofOnInflowBoundary += vertexOnInflowBoundary[dofIndex];
        }

        dirichletNodesInt[ localIndexSet.index(i)[0] ] += dofOnInflowBoundary;

      } // end dof


    } // end element e

    for(unsigned int i=0; i<dirichletNodes.size(); i++)
    {
      dirichletNodes[i] = (dirichletNodesInt[i]>0);
    }
  }

//*************************************
  template <class FEBasis, class Function>
  void BoundaryTools::getInflowBoundaryValue(
                  const FEBasis& feBasis,
                  std::vector<double>& rhsInflowContrib,
                  Function& g
                  )
{
    const unsigned int dim = FEBasis::GridView::dimension;

    typedef typename FEBasis::GridView GridView;
    GridView gridView = feBasis.gridView();

    const unsigned int dofs = feBasis.size();

    rhsInflowContrib.resize(dofs);

    auto localView = feBasis.localView();
    auto localIndexSet = feBasis.localIndexSet();

    BoundaryCondition<Function> bc = BoundaryCondition<Function>(g);
    std::vector<double> out;

    for(const auto& e : elements(gridView))
    {
      localView.bind(e);
      const auto& localFEM = localView.tree().finiteElement();
      const auto& localInterp = localFEM.localInterpolation() ;

      localIndexSet.bind(localView);

      // dofs in the current finite element
      const unsigned int dofsLocal = localFEM.localCoefficients().size();

      const unsigned int nFace
          = ReferenceElements<double, dim>::general(e.type()).size(dim-1);
      const unsigned int nVertex
          = ReferenceElements<double, dim>::general(e.type()).size(dim);

      localInterp.interpolate(bc,out);

      for(unsigned int i=0;i<dofsLocal;i++)
      {
        rhsInflowContrib[ localIndexSet.index(i)[0] ] = out[i];
      }
    }
}

//*************************************
  std::vector<unsigned int> BoundaryTools::getVertexOfIntersection(
                              unsigned int indexIntersection,
                              GeometryType geometryType
                              )
  {
    std::vector<unsigned int> indexVertex(2,0);

    if(geometryType.isSimplex()) {
      if(indexIntersection==0)
      {
        indexVertex[0]=0;
        indexVertex[1]=1;
      }
      else if(indexIntersection==1)
      {
        indexVertex[0]=0;
        indexVertex[1]=2;
      }
      else if(indexIntersection==2)
      {
        indexVertex[0]=1;
        indexVertex[1]=2;
      }
    } else if(geometryType.isCube()) {
      if(indexIntersection==0)
      {
        indexVertex[0]=0;
        indexVertex[1]=2;
      }
      else if(indexIntersection==1)
      {
        indexVertex[0]=1;
        indexVertex[1]=3;
      }
      else if(indexIntersection==2)
      {
        indexVertex[0]=0;
        indexVertex[1]=1;
      }
      else if(indexIntersection==3)
      {
        indexVertex[0]=2;
        indexVertex[1]=3;
      }
    } else {
      DUNE_THROW(Dune::NotImplemented, "getVertexOfIntersection not "
              "implemented for geometry type" << geometryType.id());
    }

    return indexVertex;

  }



} // end namespace Dune

#endif // DUNE_DPG_BOUNDARY_TOOLS
