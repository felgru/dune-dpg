#ifndef DUNE_DPG_BOUNDARY_TOOLS
#define DUNE_DPG_BOUNDARY_TOOLS

#include <iostream>

#include <array>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <dune/dpg/functions/localindexsetiteration.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dune/istl/io.hh>


namespace Dune {

  class BoundaryTools
  {

    template<class Function>
    class BoundaryCondition
    {
      const Function& g_;

    public:

      using DomainType = FieldVector<double, 2>;
      using RangeType  = FieldVector<double, 1>;

      BoundaryCondition(const Function& g) : g_(g) {};

      // Remark: this signature assumes that we have a 2D scalar problem
      void evaluate(
        const DomainType& x,
        RangeType& y) const;
    };

  public:
    BoundaryTools() = delete;

    template <class FEBasis,class Direction>
    static void getInflowBoundaryMask(
                  const FEBasis& ,
                  std::vector<bool>& ,
                  const Direction&
                  );

    template <class FEBasis>
    static void getBoundaryMask(
                  const FEBasis& ,
                  std::vector<bool>&
                  );

    template <class FEBasis, class Function>
    static void getBoundaryValue(
                  const FEBasis& ,
                  std::vector<double>& ,
                  const Function&
                  );

  private:
    static std::array<unsigned int, 2> getVerticesOfIntersection(
                        unsigned int ,
                        GeometryType
                        );
  };

  template<class Function>
  void BoundaryTools::BoundaryCondition<Function>::evaluate(
                                            const DomainType& x,
                                            RangeType& y
                                            ) const
  {
    y = g_(x);
  }

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

    const size_t dofs = feBasis.size();

    dirichletNodes.resize(dofs);
    std::vector<unsigned int> dirichletNodesInt(dofs,0);

    auto localView = feBasis.localView();
    auto localIndexSet = feBasis.localIndexSet();
    using LocalIndexSet = decltype(localIndexSet);


    for(const auto& e : elements(gridView))
    {
      localView.bind(e);
      const auto& localFE = localView.tree().finiteElement();

      localIndexSet.bind(localView);

#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
      const unsigned int nFace
          = referenceElement<double, dim>(e.type()).size(dim-1);
      const unsigned int nVertex
          = referenceElement<double, dim>(e.type()).size(dim);
#else
      const unsigned int nFace
          = ReferenceElements<double, dim>::general(e.type()).size(dim-1);
      const unsigned int nVertex
          = ReferenceElements<double, dim>::general(e.type()).size(dim);
#endif

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
          std::array<unsigned int, 2> vertexOfIntersection
              = getVerticesOfIntersection(indexIntersection, e.type());

          vertexOnInflowBoundary[ vertexOfIntersection[0] ] += 1;
          vertexOnInflowBoundary[ vertexOfIntersection[1] ] += 1;
        }
      }

      // For each dof, we check whether it belongs to the inflow boundary
      using size_type = typename LocalIndexSet::size_type;
      using MultiIndex = typename LocalIndexSet::MultiIndex;
      iterateOverLocalIndexSet(localIndexSet,
        [&](size_type i, MultiIndex gi)
        {
          unsigned int dofOnInflowBoundary = 0;

          // localkey of dof i
          const auto& dofLocalKey = localFE.localCoefficients().localKey(i);

          // Codimension and subentity index of the current dof
          const unsigned int dofCodim = dofLocalKey.codim();
          const unsigned int dofIndex = dofLocalKey.subEntity();

          if(dofCodim == 1) //the dof belongs to a face
          {
            dofOnInflowBoundary = faceOnInflowBoundary[dofIndex];
          }
          if(dofCodim == 2) //the dof belongs to a vertex
          {
            dofOnInflowBoundary = vertexOnInflowBoundary[dofIndex];
          }

          dirichletNodesInt[ gi[0] ] += dofOnInflowBoundary;

        },
        [](size_type) {},
        [&](size_type i, MultiIndex gi, double wi) {
          unsigned int dofOnInflowBoundary = 0;

          // localkey of dof i
          const auto& dofLocalKey = localFE.localCoefficients().localKey(i);

          // Codimension and subentity index of the current dof
          const unsigned int dofCodim = dofLocalKey.codim();
          const unsigned int dofIndex = dofLocalKey.subEntity();

          if(dofCodim == 1) //the dof belongs to a face
          {
            dofOnInflowBoundary = faceOnInflowBoundary[dofIndex];
          }
          if(dofCodim == 2) //the dof belongs to a vertex
          {
            dofOnInflowBoundary = vertexOnInflowBoundary[dofIndex];
          }

          dirichletNodesInt[ gi[0] ] += wi * dofOnInflowBoundary;
        });

    } // end element e

    for(size_t i=0; i<dirichletNodes.size(); i++)
    {
      dirichletNodes[i] = (dirichletNodesInt[i]>0);
    }
  }

  /**
   * \brief Writes in the vector dirichletNodes whether a degree of freedom is in the boundary (value 1) or not (value 0).
   *        The degrees of freedom are relative to the finite element basis feBasis.
   * \param feBasis        a finite element basis
   * \param dirichletNodes the vector where we store the output
   */
  template <class FEBasis>
  void BoundaryTools::getBoundaryMask(
                        const FEBasis& feBasis,
                        std::vector<bool>& dirichletNodes
                        )
  {
    const unsigned int dim = FEBasis::GridView::dimension;

    typedef typename FEBasis::GridView GridView;
    GridView gridView = feBasis.gridView();

    const size_t dofs = feBasis.size();

    dirichletNodes.resize(dofs);
    std::vector<unsigned int> dirichletNodesInt(dofs,0);

    auto localView = feBasis.localView();
    auto localIndexSet = feBasis.localIndexSet();
    using LocalIndexSet = decltype(localIndexSet);


    for(const auto& e : elements(gridView))
    {
      localView.bind(e);
      const auto& localFE = localView.tree().finiteElement();

      localIndexSet.bind(localView);

#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
      const unsigned int nFace
          = referenceElement<double, dim>(e.type()).size(dim-1);
      const unsigned int nVertex
          = referenceElement<double, dim>(e.type()).size(dim);
#else
      const unsigned int nFace
          = ReferenceElements<double, dim>::general(e.type()).size(dim-1);
      const unsigned int nVertex
          = ReferenceElements<double, dim>::general(e.type()).size(dim);
#endif

      // For every vertex, we have to see whether it is on the boundary.
      // If vertex i is on the boundary, we will have vertexOnBoundary[i] > 0.
      std::vector<unsigned int> vertexOnBoundary(nVertex, 0);

      // for all intersections, we see which one lies on the boundary
      // if intersection i lies on the boundary, then faceOnBoundary[i] == true
      // we will assume that an intersection is simply a face for us
      std::vector<unsigned int> faceOnBoundary(nFace, 0);

      for (auto&& intersection : intersections(gridView, e))
      {
        // Local index of the intersection
        const unsigned int indexIntersection = intersection.indexInInside();

        faceOnBoundary[indexIntersection] = intersection.boundary();

        // if the intersection is on the inflow boundary, we have to update
        // what are the local vertices that are also on the inflow boundary
        if(faceOnBoundary[indexIntersection])
        {
          // We see what are the vertices associated to the current
          // intersection (assumed to be a face)
          // TODO: That might give false indices on elements with hanging nodes.
          std::array<unsigned int, 2> vertexOfIntersection
              = getVerticesOfIntersection(indexIntersection, e.type());

          vertexOnBoundary[ vertexOfIntersection[0] ] += 1;
          vertexOnBoundary[ vertexOfIntersection[1] ] += 1;
        }
      }

      // For each dof, we check whether it belongs to the boundary
      using size_type = typename LocalIndexSet::size_type;
      using MultiIndex = typename LocalIndexSet::MultiIndex;
      iterateOverLocalIndexSet(localIndexSet,
        [&](size_type i, MultiIndex gi)
        {
          unsigned int dofOnBoundary = 0;

          // localkey of dof i
          const auto& dofLocalKey = localFE.localCoefficients().localKey(i);

          // Codimension and subentity index of the current dof
          const unsigned int dofCodim = dofLocalKey.codim();
          const unsigned int dofIndex = dofLocalKey.subEntity();

          if(dofCodim == 1) //the dof belongs to a face
          {
            dofOnBoundary = faceOnBoundary[dofIndex];
          }
          if(dofCodim == 2) //the dof belongs to a vertex
          {
            dofOnBoundary = vertexOnBoundary[dofIndex];
          }

          dirichletNodesInt[ gi[0] ] += dofOnBoundary;
        },
        [](size_type) {},
        [&](size_type i, MultiIndex gi, double wi) {
          unsigned int dofOnBoundary = 0;

          // localkey of dof i
          const auto& dofLocalKey = localFE.localCoefficients().localKey(i);

          // Codimension and subentity index of the current dof
          const unsigned int dofCodim = dofLocalKey.codim();
          const unsigned int dofIndex = dofLocalKey.subEntity();

          if(dofCodim == 1) //the dof belongs to a face
          {
            dofOnBoundary = faceOnBoundary[dofIndex];
          }
          if(dofCodim == 2) //the dof belongs to a vertex
          {
            dofOnBoundary = vertexOnBoundary[dofIndex];
          }

          dirichletNodesInt[ gi[0] ] += dofOnBoundary;
        });

    } // end element e

    for(size_t i=0; i<dirichletNodes.size(); i++)
    {
      dirichletNodes[i] = (dirichletNodesInt[i]>0);
    }
  }

  /**
   * \brief Interpolates a given function to a finite element space
   *
   *        This interpolated data can then be used for handling boundary
   *        values in the FE discretization.
   *
   * \todo We actually only need to interpolate at the boundary if we want
   *       to interpolate Dirichlet boundary values.
   *
   * \param feBasis     a finite element basis
   * \param rhsContrib  the vector where we store the output
   * \param g           the right hand side function to interpolate
   */
  template <class FEBasis, class Function>
  void BoundaryTools::getBoundaryValue(
                  const FEBasis& feBasis,
                  std::vector<double>& rhsContrib,
                  const Function& g
                  )
  {
    rhsContrib.resize(feBasis.size());

    auto localView = feBasis.localView();
    auto localIndexSet = feBasis.localIndexSet();
    using LocalIndexSet = decltype(localIndexSet);

    auto localG = localFunction(Functions::makeGridViewFunction(g,
                                                feBasis.gridView()));
    std::vector<double> out;

    const auto gridView = feBasis.gridView();
    for(const auto& e : elements(gridView))
    {
      localView.bind(e);
      localIndexSet.bind(localView);
      localG.bind(e);

      localView.tree().finiteElement().localInterpolation()
               .interpolate(BoundaryCondition<decltype(localG)>(localG), out);

      using size_type = typename LocalIndexSet::size_type;
      using MultiIndex = typename LocalIndexSet::MultiIndex;
      iterateOverLocalIndexSet(localIndexSet,
        [&](size_type i, MultiIndex gi) {
          rhsContrib[ gi[0] ] = out[i];
        },
        [](size_type) {},
        [&](size_type i, MultiIndex gi, double wi) {
          rhsContrib[ gi[0] ] = wi * out[i];
        });
    }
  }

  std::array<unsigned int, 2> BoundaryTools::getVerticesOfIntersection(
                              unsigned int indexIntersection,
                              GeometryType geometryType
                              )
  {
    std::array<unsigned int, 2> indexVertex{0, 0};

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
      DUNE_THROW(Dune::NotImplemented, "getVerticesOfIntersection not "
              "implemented for geometry type" << geometryType.id());
    }

    return indexVertex;

  }



} // end namespace Dune

#endif // DUNE_DPG_BOUNDARY_TOOLS
