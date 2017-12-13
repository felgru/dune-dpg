// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_QKTRANSPORTLOCALBASIS_HH
#define DUNE_QKTRANSPORTLOCALBASIS_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/power.hh>
#include <dune/common/function.hh>
#include <dune/common/version.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/localfunctions/lagrange/pk2d.hh>
#include <dune/localfunctions/lagrange/prismp1.hh>

#include <numeric>

namespace Dune
{
/**@ingroup LocalBasisImplementation
         \brief Lagrange shape functions with transport direction

         \tparam D Type to represent the field in the domain.
         \tparam R Type to represent the field in the range.
         \tparam dim Dimension of the cube

         \nosubgrouping
 */
template<class D, class R, int dim,int k>
class QkTransportLocalBasis
{
public:
    typedef LocalBasisTraits<D,dim,Dune::FieldVector<D,dim>,R,1,Dune::FieldVector<R,1>,
            Dune::FieldMatrix<R,1,dim> > Traits;

    QkTransportLocalBasis(FieldVector<D,dim> transport):beta_(transport)
    {}


    //! \brief number of shape functions
    unsigned int size () const
    {
        if(beta_[0]==0||beta_[1]==0)
            return sizeQ;
        if(beta_[0]==beta_[1])
            return sizeQ;
        if(dim==2)
            return (sizeP+sizeQ-(k+1));
        else
            static_assert(dim==2, "QkTransportLocalBasis is only implemented in 2D.");


    }

    void setCorners (std::vector<typename Traits::DomainType>& cornerT,std::vector<typename Traits::DomainType>& cornerQ,bool& triangleIsBottom) const
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

            triangleIsBottom=true;
        }
        else if(beta_[0]==0){
            cornerT={{0,0},{0,1},beta_};
            cornerQ={{0,0},{1,0},beta_,{1,1}}; //again the order is VERY IMPORTANT
            //the shape function have to be in specific order

            triangleIsBottom=false;

        }
        else{

            cornerT={{0,0},{0,1},beta_};
            cornerQ={{0,0},beta_,{1,0},{1,1}};

            triangleIsBottom=false;
        }
    }


    //! \brief Evaluate all shape functions
    // This methods figures out, where the point (that needs to be evaluated) is
    // located and then chooses the shape functions of the triagle
    //        or the rectangle
    inline void evaluateFunction (const typename Traits::DomainType& in,
            std::vector<typename Traits::RangeType>& out) const
    {
        out.resize(size());


        bool triangleIsBottom;
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
        const GeometryType triangle = GeometryTypes::triangle;
        const GeometryType quad = GeometryTypes::quadrilateral;
#else
        GeometryType triangle;
        GeometryType quad;

        triangle.makeSimplex(dim);
        quad.makeCube(dim);
#endif

        std::vector<Dune::FieldVector<D,dim>> cornerT(dim+1);
        std::vector<Dune::FieldVector<D,dim>> cornerQ(4*(dim-1));
        setCorners(cornerT,cornerQ,triangleIsBottom);

        for (size_t i=0; i<size(); i++) {

            out[i] = 0.;
        }

        if(dim==2){


            //special case: transportation direction coincide with grid
            if(beta_[0]==0||beta_[1]==0){
                std::vector<FieldVector<D,1> > valuesQ;
                Dune::MultiLinearGeometry<D,2,2> part_quad(quad,cornerQ);
                Dune::FieldVector<D,dim> plocalQ=part_quad.local(in);

                localFiniteElementQ.localBasis().evaluateFunction(plocalQ, valuesQ);

                for (size_t i=0; i<valuesQ.size(); i++) {

                    out[i]=valuesQ[i];
                }

            }
            //special case: partition leads tp two triangles
            else if(beta_[0]==beta_[1]){

                double over=(in[0]-in[1]);

                if(over>=0.){


                    std::vector<FieldVector<D,1> > valuesT;
                    Dune::MultiLinearGeometry<D,2,2> part_tria(triangle,cornerT);
                    Dune::FieldVector<D,dim> plocalT=part_tria.local(in);
                    localFiniteElementP.localBasis().evaluateFunction(plocalT, valuesT);


                    size_t i=0;
                    size_t r=0;
                    for (size_t m=0;m<k+1; m++) {
                        for (size_t j=0; j<(k+1-m); j++) {
                            out[i+r]=valuesT[i];
                            i++;
                        }
                        r+=(m+1);
                    }
                }
                else{



                    std::vector<FieldVector<D,1> > valuesT2;
                    std::vector<Dune::FieldVector<D,dim>> cornerT2(dim+1);
                    cornerT2={beta_,{0,1},{0,0}};
                    Dune::MultiLinearGeometry<D,2,2> part_tria2(triangle,cornerT2);
                    Dune::FieldVector<D,dim> plocalT2=part_tria2.local(in);
                    localFiniteElementP.localBasis().evaluateFunction(plocalT2, valuesT2);


                    size_t i=0;
                    size_t r=0;
                    for (size_t m=0;m<k+1; m++) {
                        for (size_t j=0; j<(k+1-m); j++) {
                            out[(sizeQ-1)-(i+r)]=valuesT2[i];
                            i++;
                        }
                        r+=(m+1);
                    }

                }


            }
            else{
                double over=( (beta_[1]/beta_[0])*in[0]-in[1]);

                if((0. < over && triangleIsBottom)||(0.> over && !triangleIsBottom)){


                    std::vector<FieldVector<D,1> > valuesT;
                    Dune::MultiLinearGeometry<D,2,2> part_tria(triangle,cornerT);
                    Dune::FieldVector<D,dim> plocalT=part_tria.local(in);



                    localFiniteElementP.localBasis().evaluateFunction(plocalT, valuesT);


                    for (size_t i=0; i<valuesT.size(); i++) {

                        out[i]=valuesT[i];

                    }
                }
                else{



                    std::vector<FieldVector<D,1> > valuesQ;
                    Dune::MultiLinearGeometry<D,2,2> part_quad(quad,cornerQ);
                    Dune::FieldVector<D,dim> plocalQ=part_quad.local(in);



                    localFiniteElementQ.localBasis().evaluateFunction(plocalQ, valuesQ);

                    //the first k+1 elements that are shared by the triangle and the quadrilateral
                    unsigned int n=0;
                    for (size_t i=0; i<(k+1); i++) {

                        out[n]=valuesQ[i];
                        n=n+(k+1)-i;


                    }

                    for (size_t i=(k+1); i<valuesQ.size(); i++) {

                        out[sizeP+i-(k+1)]=valuesQ[i];
                    }

                }
            }
        }

    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
            std::vector<typename Traits::JacobianType>& out) const      // return value
    {
        out.resize(size());

        bool triangleIsBottom;

#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,6)
        const GeometryType triangle = GeometryTypes::triangle;
        const GeometryType quad = GeometryTypes::quadrilateral;
#else
        GeometryType triangle;
        GeometryType quad;

        triangle.makeSimplex(dim);
        quad.makeCube(dim);
#endif

        std::vector<Dune::FieldVector<D,dim>> cornerT(dim+1);
        std::vector<Dune::FieldVector<D,dim>> cornerQ(4*(dim-1));
        setCorners(cornerT,cornerQ,triangleIsBottom);

        for (size_t i=0; i<size(); i++) {

            out[i][0] = 0.;
        }

        if(dim==2){

            if(beta_[0]==0||beta_[1]==0){
                Dune::MultiLinearGeometry<D,2,2> part_quad(quad,cornerQ);
                Dune::FieldVector<D,dim> plocalQ=part_quad.local(in);
                const auto& jacobianQ = part_quad.jacobianInverseTransposed(plocalQ);
                std::vector<FieldMatrix<D,1,dim> > referenceGradientsQ;
                localFiniteElementQ.localBasis().evaluateJacobian(plocalQ, referenceGradientsQ);
                std::vector<FieldVector<D,dim> > gradientsQ(referenceGradientsQ.size());
                for (size_t i=0; i<gradientsQ.size(); i++)
                    jacobianQ.mv(referenceGradientsQ[i][0], gradientsQ[i]);


                for (size_t i=0; i<gradientsQ.size(); i++) {
                    for (int j=0; j<dim; j++) {

                        out[i][0][j]=gradientsQ[i][j];

                    }
                }

            }
            else if(beta_[0]==beta_[1]){

                double over=(in[0]-in[1]);

                if(over>=0.){

                    Dune::MultiLinearGeometry<D,2,2> part_tria(triangle,cornerT);
                    Dune::FieldVector<D,dim> plocalT=part_tria.local(in);
                    const auto& jacobianT = part_tria.jacobianInverseTransposed(plocalT);
                    std::vector<FieldMatrix<D,1,dim> > referenceGradientsT;
                    localFiniteElementP.localBasis().evaluateJacobian(plocalT, referenceGradientsT);
                    std::vector<FieldVector<D,dim> > gradientsT(referenceGradientsT.size());
                    for (size_t i=0; i<gradientsT.size(); i++)
                        jacobianT.mv(referenceGradientsT[i][0], gradientsT[i]);


                    size_t i=0;
                    size_t r=0;
                    for (size_t m=0;m<k+1; m++) {
                        for (size_t s=0; s<(k+1-m); s++) {
                            for (int j=0; j<dim; j++) {

                                out[i+r][0][j]=gradientsT[i][j];

                            }
                            i++;
                        }
                        r+=(m+1);
                    }
                }
                else{

                    std::vector<Dune::FieldVector<D,dim>> cornerT2(dim+1);
                    cornerT2={beta_,{0,1},{0,0}};
                    Dune::MultiLinearGeometry<D,2,2> part_tria2(triangle,cornerT2);
                    Dune::FieldVector<D,dim> plocalT2=part_tria2.local(in);
                    const auto& jacobianT2 = part_tria2.jacobianInverseTransposed(plocalT2);
                    std::vector<FieldMatrix<D,1,dim> > referenceGradientsT2;
                    localFiniteElementP.localBasis().evaluateJacobian(plocalT2, referenceGradientsT2);
                    std::vector<FieldVector<D,dim> > gradientsT2(referenceGradientsT2.size());
                    for (size_t i=0; i<gradientsT2.size(); i++)
                        jacobianT2.mv(referenceGradientsT2[i][0], gradientsT2[i]);


                    size_t i=0;
                    size_t r=0;
                    for (size_t m=0;m<k+1; m++) {
                        for (size_t s=0; s<(k+1-m); s++) {
                            for (int j=0; j<dim; j++) {

                                out[(sizeQ-1)-(i+r)][0][j]=gradientsT2[i][j];

                            }
                            i++;
                        }
                        r+=(m+1);
                    }

                }


            }
            else{

                double over=( (beta_[1]/beta_[0])*in[0]-in[1]);

                if((0. < over && triangleIsBottom)||(0.> over && !triangleIsBottom)){



                    Dune::MultiLinearGeometry<D,2,2> part_tria(triangle,cornerT);
                    Dune::FieldVector<D,dim> plocalT=part_tria.local(in);
                    const auto& jacobianT = part_tria.jacobianInverseTransposed(plocalT);
                    std::vector<FieldMatrix<D,1,dim> > referenceGradientsT;
                    localFiniteElementP.localBasis().evaluateJacobian(plocalT, referenceGradientsT);
                    std::vector<FieldVector<D,dim> > gradientsT(referenceGradientsT.size());
                    for (size_t i=0; i<gradientsT.size(); i++)
                        jacobianT.mv(referenceGradientsT[i][0], gradientsT[i]);

                    for (size_t i=0; i<gradientsT.size(); i++) {
                        for (int j=0; j<dim; j++) {

                            out[i][0][j]=gradientsT[i][j];

                        }
                    }

                }
                else{

                    Dune::MultiLinearGeometry<D,2,2> part_quad(quad,cornerQ);
                    Dune::FieldVector<D,dim> plocalQ=part_quad.local(in);
                    const auto& jacobianQ = part_quad.jacobianInverseTransposed(plocalQ);
                    std::vector<FieldMatrix<D,1,dim> > referenceGradientsQ;
                    localFiniteElementQ.localBasis().evaluateJacobian(plocalQ, referenceGradientsQ);
                    std::vector<FieldVector<D,dim> > gradientsQ(referenceGradientsQ.size());
                    for (size_t i=0; i<gradientsQ.size(); i++)
                        jacobianQ.mv(referenceGradientsQ[i][0], gradientsQ[i]);


                    unsigned int n=0;
                    for (size_t i=0; i<(k+1); i++) {
                        for (int j=0; j<dim; j++) {

                            out[n][0][j]=gradientsQ[i][j];

                        }
                        n=n+(k+1)-i;

                    }

                    for (size_t i=(k+1); i<gradientsQ.size(); i++) {
                        for (int j=0; j<dim; j++) {

                            out[sizeP+i-(k+1)][0][j]=gradientsQ[i][j];

                        }
                    }


                }
            }
        }

    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    void partial(const std::array<unsigned int,d>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else {
        DUNE_THROW(Dune::NotImplemented,
            "partial only implemented for derivatives of order 0!");
      }
    }

    //! \brief Polynomial order of the shape functions
    unsigned int order () const
    {
        return k;
    }
private:
    Pk2DLocalFiniteElement<D,D,k> localFiniteElementP;
    const unsigned int sizeP= localFiniteElementP.size();
    QkLocalFiniteElement<D,D,dim,k> localFiniteElementQ;
    const unsigned int sizeQ= localFiniteElementQ.size();
    FieldVector<D,dim> beta_;
};
}
#endif
