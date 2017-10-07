// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_TRANSPORTQUADRATURERULE_HH_
#define DUNE_GEOMETRY_TRANSPORTQUADRATURERULE_HH_

#include <dune/common/version.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/type.hh>

namespace Dune {

template<typename ct, int dim>
class TransportQuadratureRule : public QuadratureRule<ct,dim> {


private:

        void setCorners (std::vector<Dune::FieldVector<ct,dim>>& cornerT,std::vector<Dune::FieldVector<ct,dim>>& cornerQ,FieldVector<ct,dim> beta_) const
        {
                if((dim==2 && beta_[1]<1) || beta_[0]==beta_[1])
                {
                        cornerT[0]={0,0};
                        cornerT[1]={1,0};
                        cornerT[2]=beta_;
                        cornerQ[0]={0,0};
                        cornerQ[1]=beta_;           //order of corners is VERY IMPORTANT
                        cornerQ[2]={0,1};                //for the output of QkFiniteElement
                        cornerQ[3]={1,1};
                }
                else if(beta_[0]==0){
                        cornerT={{0,0},{0,1},beta_};
                        cornerQ={{0,0},{1,0},beta_,{1,1}};

                }
                else{

                        cornerT={{0,0},{0,1},beta_};
                        cornerQ={{0,0},beta_,{1,0},{1,1}};
                }
        }


        void makeRule(FieldVector<ct,dim> beta) {

                std::vector<Dune::FieldVector<ct,dim>> cornerT(dim+1);
                std::vector<Dune::FieldVector<ct,dim>> cornerQ(4*(dim-1));
                setCorners(cornerT,cornerQ,beta);

#if DUNE_VERSION_NEWER(DUNE_GRID,2,6)
                MultiLinearGeometry<ct,2,2> part_tria(GeometryTypes::triangle, cornerT);
#else
                GeometryType triangle;
                triangle.makeSimplex(dim);
                MultiLinearGeometry<ct,2,2> part_tria(triangle,cornerT);
#endif
                const QuadratureRule<ct,dim>& quadP = Dune::template QuadratureRules<ct,dim>::rule(part_tria.type(), this->delivered_order);
                for (unsigned int pt=0; pt < quadP.size(); ++pt) {// loop over quadrature points

                        const Dune::FieldVector<ct,dim>& quadPos = quadP[pt].position(); // local position of the quadrature point in the refined element
                        const Dune::FieldVector<ct,dim>& global = part_tria.global(quadPos); // global position of the qpt == local position in element


                        ct ie = (part_tria.integrationElement(quadPos) /1.);
                        ct weight = quadP[pt].weight()*ie; // weight of the quadrature point

                        this->push_back(QuadraturePoint<ct,dim>(global,weight));
                }



#if DUNE_VERSION_NEWER(DUNE_GRID,2,6)
                MultiLinearGeometry<ct,2,2> part_quad(GeometryTypes::quadrilateral, cornerQ);
#else
                GeometryType quad;
                quad.makeCube(dim);
                MultiLinearGeometry<ct,2,2> part_quad(quad,cornerQ);
#endif
                const QuadratureRule<ct,dim>& quadQ = Dune::template QuadratureRules<ct,dim>::rule(part_quad.type(), this->delivered_order);
                for (unsigned int pt=0; pt < quadQ.size(); ++pt) {// loop over quadrature points

                        const Dune::FieldVector<ct,dim>& quadPos = quadQ[pt].position(); // local position of the quadrature point in the refined element
                        const Dune::FieldVector<ct,dim>& global = part_quad.global(quadPos); // global position of the qpt == local position in element


                        ct ie = (part_quad.integrationElement(quadPos) /1.);
                        ct weight = quadQ[pt].weight()*ie; // weight of the quadrature point

                        this->push_back(QuadraturePoint<ct,dim>(global,weight));
                }


        }

public:
        TransportQuadratureRule(GeometryType t, int order,FieldVector<ct,dim> beta)
            : QuadratureRule<ct,dim>(t, order)
        {
                assert(this->geometry_type.isCube());
                makeRule(beta);
        }
};
}



#endif /* DUNE_GEOMETRY_TRANSPORTQUADRATURERULE_HH_ */
