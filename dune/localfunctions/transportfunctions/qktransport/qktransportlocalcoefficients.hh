/*
 * qktransportlocalcoefficients.hh
 *
 *  Created on: April 27, 2015
 *      Author: koenig
 */

#ifndef DUNE_QKTRANSPORTLOCALCOEFFICIENTS_HH_
#define DUNE_QKTRANSPORTLOCALCOEFFICIENTS_HH_

#include <cstddef>
#include <iostream>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/localfunctions/lagrange/pk2d.hh>

namespace Dune
{

/**@ingroup LocalLayoutImplementation
         \brief Layout map for Q1 elements

         \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
 */
template <int d,int k>
class QkTransportLocalCoefficients
{




    void setup2d(std::vector<unsigned int>& subEntity)
    {
        assert(k>0);
        unsigned lastIndex=0;

        // LocalKey: entity number, entity codim, dof indices within each entity
        /* edge and vertex numbering
                 2----3----3
                 |         |
                 |         |
                 0         1
                 |         |
                 |         |
                 0----2----1
         */
        if(beta_[0]==1.){                                                //triangle at bottom
            // lower edge (2)
            subEntity[lastIndex++] = 0;                 // corner 0
            for (unsigned i = 0; i < k - 1; ++i)
                subEntity[lastIndex++] = 2;           // inner dofs of lower edge (2)

            subEntity[lastIndex++] = 1;                 // corner 1


            for (unsigned int j=0; j<=k-1; j++)              //for the triangle part
                for (unsigned int i=0; i<=k-1-j; i++)
                {
                    if (i+j==(k-1))
                        subEntity[lastIndex++] = 1;
                    else
                        subEntity[lastIndex++] = 0;
                }

            // iterate from bottom to top over inner edge dofs
            for (unsigned e = 0; e < k - 1; ++e) {
                subEntity[lastIndex++] = 0;                   // left edge (0)
                for (unsigned i = 0; i < k - 1; ++i)
                    subEntity[lastIndex++] = 0;                     // face dofs
                subEntity[lastIndex++] = 1;                   // right edge (1)
            }

            // upper edge (3)
            subEntity[lastIndex++] = 2;                 // corner 2
            for (unsigned i = 0; i < k - 1; ++i)
                subEntity[lastIndex++] = 3;                   // inner dofs of upper edge (3)

            subEntity[lastIndex++] = 3;                 // corner 3
        }
        else{ //triangle up
            // left edge (0)
            subEntity[lastIndex++] = 0;                 // corner 0
            for (unsigned i = 0; i < k - 1; ++i)
                subEntity[lastIndex++] = 0;           // inner dofs of lower edge (0)

            subEntity[lastIndex++] = 3;                 // corner 3


            for (unsigned int j=0; j<=k-1; j++)              //for the triangle part
                for (unsigned int i=0; i<=k-1-j; i++)
                {
                    if (i+j==(k-1))
                        subEntity[lastIndex++] = 3;
                    else
                        subEntity[lastIndex++] = 0;
                }

            // iterate from left to right over inner edge dofs
            for (unsigned e = 0; e < k - 1; ++e) {
                subEntity[lastIndex++] = 2;                   // lower edge (2)
                for (unsigned i = 0; i < k - 1; ++i)
                    subEntity[lastIndex++] = 0;                     // face dofs
                subEntity[lastIndex++] = 3;                   // upper edge (3)
            }

            // upper edge (3)
            subEntity[lastIndex++] = 1;                 // corner 1
            for (unsigned i = 0; i < k - 1; ++i)
                subEntity[lastIndex++] = 1;                   // inner dofs of right edge (1)

            subEntity[lastIndex++] = 3;                 // corner 3

        }

        assert(((sizeP+sizeQ-(k+1))==lastIndex));


    }




public:
    QkTransportLocalCoefficients (FieldVector<double,d> transport)
    {
        beta_=transport;

        li.resize(size());



        auto coefficientP=localFiniteElementP.localCoefficients();
        auto coefficientQ=localFiniteElementQ.localCoefficients();


        if(beta_[0]==0||beta_[1]==0||beta_[0]==beta_[1])
            for (size_t i=0; i<size(); i++)
                li[i] = coefficientQ.localKey(i);
        else{

            // Set up array of codimension-per-dof-number
            std::vector<unsigned int> codim(size());

            unsigned int n=0;
            for (std::size_t i=0; i<sizeP; i++) {
                codim[i] = 0;

                int codimInPk= coefficientP.localKey(i).codim();

                codim[i]=codimInPk;
            }
            for (std::size_t i=0; i<k-1; i++) {
                n=n+(k+1)-i;

                codim[n]=0;
            }
            codim[sizeP-1]=1;

            for (std::size_t i=sizeP; i<codim.size(); i++) {
                codim[i] = 0;

                int codimInQk= coefficientQ.localKey((i-sizeP+(k+1))).codim();
                codim[i]=codimInQk;
            }

            // Set up array of index for each subEntity
            std::vector<unsigned int> index(size());
            unsigned int c=0;
            n=0;
            for (std::size_t i=0; i<k; i++) {
                index[i] = 0;

                int indexInPk= coefficientP.localKey(i).index();

                index[i]=indexInPk;
            }
            for (std::size_t i=k+1; i<sizeP-1; i++) {

                if(coefficientP.localKey(i).subEntity()==1 && coefficientP.localKey(i).codim()==1)
                {
                    index[i]=n;
                    n++;
                }
                else{
                    index[i]=c;
                    c++;
                }
            }
            index[sizeP-1]=k-1;

            for (std::size_t i=sizeP; i<index.size(); i++) {
                index[i] = 0;
                auto currentLocalKeyQ =coefficientQ.localKey(i-sizeP+(k+1));
                index[i]=currentLocalKeyQ.index();
                if(currentLocalKeyQ.subEntity()==0 && currentLocalKeyQ.codim()==0)
                    index[i]=(currentLocalKeyQ.index()+(k*(k-1)/2));
                else if(currentLocalKeyQ.subEntity()==1 && currentLocalKeyQ.codim()==1)
                    index[i]=(currentLocalKeyQ.index()+k);
                else
                    index[i]=currentLocalKeyQ.index();
            }


            // Set up entity and dof numbers for each (supported) dimension separately
                    std::vector<unsigned int> subEntity(size());

            if (k==1) {

                if(d==2 && beta_[0]==1.){
                    for (std::size_t i=0; i<sizeP; i++)
                        subEntity[i] = i;
                    subEntity[2]=1;
                    for (std::size_t i=sizeP; i<size(); i++)
                        subEntity[i] = i-1;
                }
                else{
                    subEntity[0] = 0;
                    subEntity[1] = 2;
                    subEntity[2] = 3;
                    subEntity[3] = 1;
                    subEntity[4] = 3;
                }


            } else if (d==2) {

                setup2d(subEntity);

                /*} else if (d==3) {

                setup3d(subEntity);*/

                } else
                    DUNE_THROW(Dune::NotImplemented, "QkLocalCoefficients for k==" << k << " and d==" << d);

            for (size_t i=0; i<size(); i++)
                li[i] = LocalKey(subEntity[i], codim[i], index[i]);
        }
    }

    //! number of coefficients
    std::size_t size () const
    {
        if(beta_[0]==0||beta_[1]==0||beta_[0]==beta_[1])
            return sizeQ;
        else
            return (sizeP+sizeQ-(k+1));
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
        return li[i];
    }

private:
    std::vector<LocalKey> li;
    Pk2DLocalFiniteElement<double,double,k> localFiniteElementP;
    const unsigned int sizeP= localFiniteElementP.size();
    QkLocalFiniteElement<double,double,d,k> localFiniteElementQ;
    const unsigned int sizeQ= localFiniteElementQ.size();
    FieldVector<double,d> beta_;

};

}




#endif /* DUNE_QKTRANSPORTLOCALCOEFFICIENTS_HH_ */
