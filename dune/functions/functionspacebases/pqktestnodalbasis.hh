/*
 * pqktestnodalbasis.hh
 *
 *  Created on: April 29, 2015
 *      Author: koenig
 */

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ1TESTBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ1TestBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>
#include <dune/common/std/final.hh>

#include <dune/localfunctions/transportfunctions/pqktestfactory.hh>
#include <dune/localfunctions/transportfunctions/qktest.hh>

#include <dune/typetree/leafnode.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/localfunctions/lagrange/pk2d.hh>

#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>


namespace Dune {
namespace Functions {

template<typename GV, int k>
class PQkTestBasisLocalView;

template<typename GV, int k>
class PQkTestBasisLeafNode;

template<typename GV, int k>
class PQkTestIndexSet;

template<typename GV, int k>
class PQkTestLocalIndexSet
{
        enum {dim = GV::dimension};

public:
        typedef std::size_t size_type;

        /** \brief Type of the local view on the restriction of the basis to a single
element */
        typedef PQkTestBasisLocalView<GV,k> LocalView;

        /** \brief Type used for global numbering of the basis vectors */
        typedef std::array<size_type, 1> MultiIndex;

        PQkTestLocalIndexSet(const PQkTestIndexSet<GV,k> & indexSet)
        : basisIndexSet_(indexSet)
        {}

        /** \brief Bind the view to a grid element
         *
         * Having to bind the view to an element before being able to actually access any
of its data members
         * offers to centralize some expensive setup code in the 'bind' method, which can
save a lot of run-time.
         */
        void bind(const PQkTestBasisLocalView<GV,k>& localView)
        {
                localView_ = &localView;
        }

        /** \brief Unbind the view
         */
        void unbind()
        {
                localView_ = nullptr;
        }

        /** \brief Size of subtree rooted in this node (element-local)
         */
        size_type size() const
        {

                return localView_->tree().finiteElement_->size();

        }

        //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
        const MultiIndex index(size_type i) const
        {
                const auto& gridIndexSet = basisIndexSet_.gridView_.indexSet();
                const auto& element = localView_->element();

                const auto& beta = basisIndexSet_.transport();


                if(beta[0]==0||beta[1]==0||beta[0]==beta[1])
                    return {basisIndexSet_.sizeQ*gridIndexSet.subIndex(element,0,0) + i};
                else
                    return {basisIndexSet_.dofsPerQuad*gridIndexSet.subIndex(element,0,0) + i};
                //DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");

                //DUNE_THROW(Dune::NotImplemented, "No index method for " << dim << "d grids available yet!");
        }

        /** \brief Return the local view that we are attached to
         */
        const LocalView& localView() const
        {
                return *localView_;
        }

        const PQkTestBasisLocalView<GV,k>* localView_;

        const PQkTestIndexSet<GV,k> basisIndexSet_;


};

template<typename GV, int k>
class PQkTestIndexSet
{
        static const int dim = GV::dimension;

        // Needs the mapper
        friend class PQkTestLocalIndexSet<GV,k>;


public:

        typedef PQkTestLocalIndexSet<GV,k> LocalIndexSet;

        PQkTestIndexSet(const GV& gridView,FieldVector<double,dim> transport)
        : gridView_(gridView),
          beta_(transport)
        {}

        std::size_t size() const
        {

                GeometryType quad;
                quad.makeQuadrilateral();
                if(beta_[0]==0||beta_[1]==0||beta_[0]==beta_[1])
                    return sizeQ * gridView_.size(quad);
                else
                    return dofsPerQuad * gridView_.size(quad);
                //DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
        }

        LocalIndexSet localIndexSet() const
        {
                return LocalIndexSet(*this);
        }

        FieldVector<double,dim> transport() const
        {
                       return beta_;
        }

private:

        //size_t quadrilateralOffset_;

        const GV gridView_;
        // Precompute the number of dofs per entity type
        Pk2DLocalFiniteElement<double,double,k> localFiniteElementP;
        const unsigned int sizeP= localFiniteElementP.size();
        QkLocalFiniteElement<double,double,dim,k> localFiniteElementQ;
        const unsigned int sizeQ= localFiniteElementQ.size();
        const int dofsPerQuad        = (sizeP+sizeQ-(k+1));
        FieldVector<double,dim> beta_;
};

/** \brief Nodal basis of a scalar third-order Lagrangean finite element space
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - Grids must be 1d or 2d
 *
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k>
class PQkTestBasis
                : public GridViewFunctionSpaceBasis<GV,
                  PQkTestBasisLocalView<GV,k>,
                  PQkTestIndexSet<GV,k>,
                  std::array<std::size_t, 1>>
                  {
        static const int dim = GV::dimension;

public:

        /** \brief The grid view that the FE space is defined on */
        typedef GV GridView;
        typedef std::size_t size_type;

        /** \brief Type of the local view on the restriction of the basis to a single
element */
        typedef PQkTestBasisLocalView<GV,k> LocalView;

        /** \brief Type used for global numbering of the basis vectors */
        typedef std::array<size_type, 1> MultiIndex;

        /** \brief Constructor for a given grid view object */
        PQkTestBasis(const GridView& gv, FieldVector<double,dim> transport) :
                gridView_(gv),
                indexSet_(gv,transport),
                beta_(transport)
        {}

        /** \brief Obtain the grid view that the basis is defined on
         */
        const GridView& gridView() const DUNE_FINAL
                        {
                return gridView_;
                        }

        PQkTestIndexSet<GV,k> indexSet() const
                                                  {
                return indexSet_;
                                                  }

        /** \brief Return local view for basis
         *
         */
        LocalView localView() const
        {
                return LocalView(this);
        }

        FieldVector<double,dim> transport() const
                                                {
                return beta_;
                                                }

protected:
        const GridView gridView_;

        PQkTestIndexSet<GV,k> indexSet_;

        FieldVector<double,dim> beta_;
};


/** \brief The restriction of a finite element basis to a single element */
template<typename GV, int k>
class PQkTestBasisLocalView
{
public:
        /** \brief The global FE basis that this is a view on */
        typedef PQkTestBasis<GV,k> GlobalBasis;
        typedef typename GlobalBasis::GridView GridView;

        /** \brief The type used for sizes */
        typedef typename GlobalBasis::size_type size_type;

        /** \brief Type used to number the degrees of freedom
         *
         * In the case of mixed finite elements this really can be a multi-index, but for
a standard
         * P3 space this is only a single-digit multi-index, i.e., it is an integer.
         */
        typedef typename GlobalBasis::MultiIndex MultiIndex;

        /** \brief Type of the grid element we are bound to */
        typedef typename GridView::template Codim<0>::Entity Element;

        /** \brief Tree of local finite elements / local shape function sets
         *
         * In the case of a P3 space this tree consists of a single leaf only,
         * i.e., Tree is basically the type of the LocalFiniteElement
         */
        typedef PQkTestBasisLeafNode<GV,k> Tree;

        /** \brief Construct local view for a given global finite element basis */
        PQkTestBasisLocalView(const GlobalBasis* globalBasis) :
                globalBasis_(globalBasis),
                tree_(globalBasis)
        {}

        /** \brief Bind the view to a grid element
         *
         * Having to bind the view to an element before being able to actually access any
of its data members
         * offers to centralize some expensive setup code in the 'bind' method, which can
save a lot of run-time.
         */
        void bind(const Element& e)
        {
                element_ = &e;
                tree_.bind(e);
        }

        /** \brief Return the grid element that the view is bound to
         *
         * \throws Dune::Exception if the view is not bound to anything
         */
        const Element& element() const
        {
                if (element_)
                        return *element_;
                else
                        DUNE_THROW(Dune::Exception, "Can't query element of unbound local view");
        }

        /** \brief Unbind from the current element
         *
         * Calling this method should only be a hint that the view can be unbound.
         * And indeed, in the PQ1TestBasisView implementation this method does nothing.
         */
        void unbind()
        {}

        /** \brief Return the local ansatz tree associated to the bound entity
         *
         * \returns Tree // This is tree
         */
        const Tree& tree() const
        {
                return tree_;
        }

        /** \brief Number of degrees of freedom on this element
         */
        size_type size() const
        {
                // We have subTreeSize==lfe.size() because we're in a leaf node.
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
                return tree_.finiteElement_->size();
#else
                return tree_.finiteElement_->localBasis().size();
#endif
        }

        /**
         * \brief Maximum local size for any element on the GridView
         *
         * This is the maximal size needed for local matrices
         * and local vectors, i.e., the result is
         *
         * The method returns (k+1)^dim, which is the number of degrees of freedom you get
for cubes.
         */
        size_type maxSize() const
        {
                return StaticPower<(k+1),GV::dimension>::power;
        }

        /** \brief Return the global basis that we are a view on
         */
        const GlobalBasis& globalBasis() const
        {
                return *globalBasis_;
        }

protected:
        const GlobalBasis* globalBasis_;
        const Element* element_;
        Tree tree_;
};


template<typename GV, int k>
class PQkTestBasisLeafNode :
                public GridFunctionSpaceBasisLeafNodeInterface<
                typename GV::template Codim<0>::Entity,
                typename Dune::PQkTestLocalFiniteElementCache<typename GV::ctype, double,
                GV::dimension,k>::FiniteElementType,
                typename PQkTestBasis<GV,k>::size_type>
{
        typedef PQkTestBasis<GV,k> GlobalBasis;
        static const int dim = GV::dimension;

        typedef typename GV::template Codim<0>::Entity E;
        typedef typename Dune::PQkTestLocalFiniteElementCache<typename GV::ctype, double,
                        dim,k> FiniteElementCache;
        typedef typename FiniteElementCache::FiniteElementType FE;
        typedef typename GlobalBasis::size_type ST;
        typedef typename GlobalBasis::MultiIndex MI;

        typedef typename GlobalBasis::LocalView LocalView;

        friend LocalView;
        friend class PQkTestLocalIndexSet<GV,k>;

public:
        typedef GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST> Interface;
        typedef typename Interface::size_type size_type;
        typedef typename Interface::Element Element;
        typedef typename Interface::FiniteElement FiniteElement;

        PQkTestBasisLeafNode(const GlobalBasis* globalBasis) :
                globalBasis_(globalBasis),
                finiteElement_(nullptr),
                element_(nullptr)
        {}

        //! Return current element, throw if unbound
        const Element& element() const DUNE_FINAL
                        {
                return *element_;
                        }

        /* * \brief Return the LocalFiniteElement for the element we are bound to
         *
         * The LocalFiniteElement implements the corresponding interfaces of the
dune-localfunctions module
         */
        const FiniteElement& finiteElement() const DUNE_FINAL
                        {
                return *finiteElement_;
                        }

        /* * \brief Size of subtree rooted in this node (element-local)
         */
        size_type size() const DUNE_FINAL
                        {
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
                return finiteElement_->size();
#else
                return finiteElement_->localBasis().size();
#endif
                        }

        //! Maps from subtree index set [0..subTreeSize-1] into root index set (element-local) [0..localSize-1]
        size_type localIndex(size_type i) const DUNE_FINAL
                        {
                return i;
                        }

        void setLocalIndex(size_type leafindex, size_type localindex) DUNE_FINAL
                        {
                DUNE_THROW(NotImplemented,"not implemented");
                        }

protected:

        //! Bind to element.
        void bind(const Element& e)
        {
                element_ = &e;
                FieldVector<typename GV::ctype,dim> normalized,unnormalized;
                unnormalized=globalBasis_->transport();
                if(unnormalized[0]!=0 && unnormalized[1]!=0){

                        double factor=unnormalized[1]/unnormalized[0];
                        if(factor<=1){
                                normalized[0]=1.;
                                normalized[1]=factor;
                        }
                        else if(factor>1){
                                normalized[0]=1./factor;
                                normalized[1]=1.;
                        }

                }
                if(unnormalized[0]==0){
                        normalized[0]=0.;
                        normalized[1]=1.;
                }
                if(unnormalized[1]==0){
                        normalized[0]=1.;
                        normalized[1]=0.;
                }
                finiteElement_ = &(cache_.get(element_->type(),normalized));
        }

        const GlobalBasis* globalBasis_;
        FiniteElementCache cache_;
        const FiniteElement* finiteElement_;
        const Element* element_;
};



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQ1TESTBASIS_HH
