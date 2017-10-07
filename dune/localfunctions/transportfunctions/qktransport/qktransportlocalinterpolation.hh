/*
 * qktransportlocalinterpolation.hh
 *
 *  Created on: April 27, 2015
 *      Author: koenig
 */

#ifndef DUNE_QKTRANSPORTLOCALINTERPOLATION_HH_
#define DUNE_QKTRANSPORTLOCALINTERPOLATION_HH_

#include <vector>
#include <dune/common/fvector.hh>
#include <dune/common/power.hh>
#include <dune/common/version.hh>

#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/localfunctions/lagrange/pk2d.hh>

namespace Dune
{

/** \todo Please doc me! */
template<int dim, class LB,int k>
class QkTransportLocalInterpolation
{
public:
    typedef typename LB::Traits::DomainFieldType D;


    QkTransportLocalInterpolation(FieldVector<D,dim> transport)
    : beta_(transport)
    {}



    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
        typename LB::Traits::DomainType x;
        typename LB::Traits::DomainType local;
        typename LB::Traits::RangeType y;




        out.resize((sizeP+sizeQ-(k+1)));

#if DUNE_VERSION_NEWER(DUNE_GRID,2,6)
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
        setCorners(cornerT,cornerQ);

        std::vector<FieldVector<D,1> > valuesT;
        Dune::MultiLinearGeometry<D,2,2> part_tria(triangle,cornerT);

        std::vector<FieldVector<D,1> > valuesQ;
        Dune::MultiLinearGeometry<D,2,2> part_quad(quad,cornerQ);

        if(beta_[0]==0||beta_[1]==0||beta_[0]==beta_[1])
        {
            out.resize(sizeQ);

            for (int i=0; i< sizeQ; i++) {


                // convert index i to multiindex
                Dune::FieldVector<int,dim> alpha(multiindex(i));

                // Generate coordinate of the i-th Lagrange point
                for (int j=0; j<dim; j++)
                    local[j] = (1.0*alpha[j])/k;

                x=part_quad.global(local);
                f.evaluate(x,y); out[i] = y;

            }


        }
        else{
            for (int i=0; i< sizeP; i++) {

                for (int j=0; j<=k; j++)
                    for (int m=0; m<=k-j; m++)
                    {
                        local[0] = ((D)m)/((D)kdiv); local[1] = ((D)j)/((D)kdiv);
                        x=part_tria.global(local);
                        f.evaluate(x,y); out[i] = y;
                    }
            }
            for (int i=k+1; i< sizeQ; i++) {


                // convert index i to multiindex
                Dune::FieldVector<int,dim> alpha(multiindex(i));

                // Generate coordinate of the i-th Lagrange point
                for (int j=0; j<dim; j++)
                    local[j] = (1.0*alpha[j])/k;

                x=part_quad.global(local);
                f.evaluate(x,y); out[i] = y;

            }
        }
    }


private:

    static const int kdiv = (k == 0 ? 1 : k);

    // Return i as a d-digit number in the (k+1)-nary system
    static Dune::FieldVector<int,dim> multiindex (int i)
    {
        Dune::FieldVector<int,dim> alpha;
        for (int j=0; j<dim; j++)
        {
            alpha[j] = i % (k+1);
            i = i/(k+1);
        }
        return alpha;
    }

    void setCorners (std::vector<FieldVector<D,dim>>& cornerT,std::vector<FieldVector<D,dim>>& cornerQ) const
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

    Pk2DLocalFiniteElement<D,D,k> localFiniteElementP;
    const unsigned int sizeP= localFiniteElementP.size();
    QkLocalFiniteElement<D,D,dim,k> localFiniteElementQ;
    const unsigned int sizeQ= localFiniteElementQ.size();
    FieldVector<D,dim> beta_;

};
}




#endif /* DUNE_QKTRANSPORTLOCALINTERPOLATION_HH_ */
