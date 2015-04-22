// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_OPTIMALTESTBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_OPTIMALTESTBASIS_HH







#include <tuple>
#include <functional>
#include <memory>
#include <type_traits>

/* 7 would be enought, but better take some more, in case we
 * change our procedures later. */
#define BOOST_FUSION_INVOKE_PROCEDURE_MAX_ARITY 10

#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/adapted/array.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/container/set/convert.hpp>
#include <boost/fusion/algorithm/auxiliary/copy.hpp>
#include <boost/fusion/algorithm/transformation/join.hpp>
#include <boost/fusion/algorithm/transformation/transform.hpp>
#include <boost/fusion/algorithm/transformation/zip.hpp>
#include <boost/fusion/algorithm/iteration/accumulate.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/functional/generation/make_fused_procedure.hpp>

#include <boost/fusion/sequence/intrinsic/size.hpp>
#include <boost/fusion/algorithm/transformation/pop_back.hpp>






#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>
#include <dune/common/std/final.hh>

#include <dune/localfunctions/optimaltestfunctions/optimaltest.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>
#include <dune/dpg/assemble_helper.hh>


namespace Dune {
namespace Functions {

template<typename BilinearForm, typename InnerProduct>
class OptimalTestBasisLocalView;

template<typename BilinearForm, typename InnerProduct>
class OptimalTestBasisLeafNode;

template<typename BilinearForm, typename InnerProduct>
class OptimalTestIndexSet;

/////////////////////////////////////////////////////////////////////////////////
// Helper begin
/////////////////////////////////////////////////////////////////////////////////


struct computeIndex
{
    computeIndex(size_t& space_index, size_t& index_result, bool& index_found)
    : space_index(&space_index),
    index_result(&index_result),
    index_found(&index_found)
    {};

    template<class LIS>
    void operator()(LIS& localIndexSet) const
    {
      if (!(*index_found))
      {
        if (localIndexSet->size()>*index_result)
        {
          *index_found=true;
          *index_result=(localIndexSet->index(*index_result))[0];
        }
        else
        {
          *space_index+=1;
          *index_result-=localIndexSet->size();
        }
      }
    }

private:
    size_t* space_index;
    size_t* index_result;
    bool* index_found;
};



template<typename MatrixType>
class Cholesky
{
public:
  Cholesky<MatrixType>(MatrixType& matrix):
  matrix(matrix)
  {
    unsigned int n = matrix.N();
    if(!(n==matrix.M()))
    {
      DUNE_THROW(Dune::Exception, "Cholesky only possible for square-matrices");
    }
    for (unsigned int k=0; k<n; k++)
    {
      double summe = matrix[k][k];
      for (unsigned int i=0; i<k; ++i)
      {
        summe -= (matrix[k][i]*matrix[k][i]*matrix[i][i]);
      }
      if (summe<1e-20)
      {
        DUNE_THROW(Dune::Exception, "Matrix for Cholesky not positive semidefinite");
      }
      matrix[k][k]=summe;
      for (unsigned int j=k+1; j<n; ++j)
      {
        summe = matrix[j][k];
        for (unsigned int i=0; i<k; ++i)
        {
          summe -= (matrix[j][i]*matrix[k][i]*matrix[i][i]); //TODO matrix[k][i]*matrix[i][i] zwischenspeichern guenstiger?
        }
        matrix[j][k] = (summe/matrix[k][k]);
      }
    }
  }

  void apply(MatrixType& rhsMatrix)
  {
    unsigned int m = rhsMatrix.M();
    unsigned int n = matrix.N();
    for (unsigned int im=0; im<m; im++)
    {
    // solve LY = B
      for (int j=0; j<n; j++)
      {
        double summe = 0;
        for (unsigned int i=0; i<j; i++)
        {
          summe+=matrix[j][i]*rhsMatrix[i][im];
        }
        rhsMatrix[j][im] = (rhsMatrix[j][im]-summe);
      }
     // solve L^TX = D^(-1)B
      for (int j=n-1; j>-1; j--)
      {
        rhsMatrix[j][im] = rhsMatrix[j][im]/matrix[j][j];
        double summe = 0;
        for (unsigned int i=j+1; i<n; i++)
        {
          summe+=matrix[i][j]*rhsMatrix[i][im];
        }
        rhsMatrix[j][im] = (rhsMatrix[j][im]-summe);
      }
    }
  }

private:
  MatrixType& matrix;
};

/////////////////////////////////////////////////////////////////////////////////
// Helper end
/////////////////////////////////////////////////////////////////////////////////




template<typename BilinearForm, typename InnerProduct>
class OptimalTestLocalIndexSet
{
  typedef typename std::tuple_element<0,typename BilinearForm::SolutionSpaces>::type::GridView GridView;
  enum {dim = GridView::dimension};

public:
  typedef typename BilinearForm::SolutionSpaces SolutionSpaces;

  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef OptimalTestBasisLocalView<BilinearForm, InnerProduct> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type, 1> MultiIndex;

  OptimalTestLocalIndexSet(const OptimalTestIndexSet<BilinearForm, InnerProduct> & indexSet, const SolutionSpaces& solSp)
  : basisIndexSet_(indexSet),
  solutionBasisIndexSet_(boost::fusion::as_vector(boost::fusion::transform(solSp,getIndexSet()))),
  solutionLocalIndexSet_(boost::fusion::as_vector(transform(solutionBasisIndexSet_,getLocalIndexSet())))
  {}

  typedef typename boost::fusion::result_of::as_vector<
             typename boost::fusion::result_of::transform<SolutionSpaces,
                                                          getIndexSet>::type
             >::type SolutionBasisIndexSet;
  typedef typename boost::fusion::result_of::as_vector<
             typename boost::fusion::result_of::transform<SolutionBasisIndexSet,
                                                          getLocalIndexSet>::type
             >::type SolutionLocalIndexSet;

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const OptimalTestBasisLocalView<BilinearForm, InnerProduct>& localView)
  {
    localView_ = &localView;

    boost::fusion::for_each(boost::fusion::zip(solutionLocalIndexSet_, localView_->tree().localViewSolution),
             boost::fusion::make_fused_procedure(bindLocalIndexSet()));
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    localView_ = nullptr;
    for_each(solutionLocalIndexSet_, applyUnbind());
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
    /* set up global offsets*/
    size_t globalOffsets[std::tuple_size<SolutionSpaces>::value];
    boost::fusion::fold(boost::fusion::zip(globalOffsets, solutionBasisIndexSet_), 0, globalOffsetHelper());

    size_t space_index=0;
    size_t index_result=i;
    bool index_found=false;

    boost::fusion::for_each(solutionLocalIndexSet_, computeIndex(space_index, index_result, index_found));

    MultiIndex result;
    result[0]=(globalOffsets[space_index]+index_result);
    return result;
    DUNE_THROW(Dune::NotImplemented, "TEXT");
  }

  /** \brief Return the local view that we are attached to
   */
  const LocalView& localView() const
  {
    return *localView_;
  }

  const OptimalTestBasisLocalView<BilinearForm, InnerProduct>* localView_;

  const OptimalTestIndexSet<BilinearForm, InnerProduct> basisIndexSet_;

  SolutionBasisIndexSet solutionBasisIndexSet_;
  SolutionLocalIndexSet solutionLocalIndexSet_;
};



template<typename BilinearForm, typename InnerProduct>
class OptimalTestIndexSet
{
  typedef typename std::tuple_element<0,typename BilinearForm::SolutionSpaces>::type::GridView GridView;
  static const int dim = GridView::dimension;

  // Needs the mapper
  friend class OptimalTestLocalIndexSet<BilinearForm, InnerProduct>;

public:
  typedef typename BilinearForm::SolutionSpaces SolutionSpaces;

  typedef OptimalTestLocalIndexSet<BilinearForm, InnerProduct> LocalIndexSet;

  OptimalTestIndexSet(const SolutionSpaces& solSp)
  : solutionSpaces(solSp),
    gridView_(std::get<0>(solSp).gridView())
  {}

  std::size_t size() const
  {
    return boost::fusion::fold(boost::fusion::transform(solutionSpaces,
                                                        getIndexSetSize()),
                               0, std::plus<std::size_t>());
  }

  LocalIndexSet localIndexSet() const
  {
    return LocalIndexSet(*this, solutionSpaces);
  }

private:
  const SolutionSpaces& solutionSpaces;
  const GridView gridView_;
};



/** \brief Nodal basis of a scalar third-order Lagrangean finite element space
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - Grids must be 1d, 2d, or 3d
 * - 3d grids must be simplex grids
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */

template<typename BilinearForm, typename InnerProduct>
class OptimalTestBasis
: public GridViewFunctionSpaceBasis<typename std::tuple_element<0,typename BilinearForm::SolutionSpaces>
                                                  ::type::GridView,
                                    OptimalTestBasisLocalView<BilinearForm,InnerProduct>,
                                    OptimalTestIndexSet<BilinearForm,InnerProduct>,
                                    std::array<std::size_t, 1> >
{
public:
  typedef typename std::tuple_element<0,typename BilinearForm::SolutionSpaces>::type::GridView GridView;

private:
  static const int dim = GridView::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef OptimalTestBasisLocalView<BilinearForm,InnerProduct> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type, 1> MultiIndex;   //TODO: Brauche ich das?

  typedef typename BilinearForm::SolutionSpaces SolutionSpaces;
  typedef typename BilinearForm::TestSpaces EnrichedTestspaces;

  /** \brief Constructor for a given grid view object */
  OptimalTestBasis(BilinearForm& bilinForm, InnerProduct& inProd) :
    gridView_(std::get<0>(bilinForm.getSolutionSpaces()).gridView()),
    bilinearForm(bilinForm),
    innerProduct(inProd),
    indexSet_(bilinearForm.getSolutionSpaces())
  {}

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const DUNE_FINAL
  {
    return gridView_;
  }

  OptimalTestIndexSet<BilinearForm, InnerProduct> indexSet() const
  {
    return indexSet_;
  }

  /** \brief Return local view for basis
   *
   */
  LocalView localView() const
  {
    return LocalView(this, bilinearForm, innerProduct);
  }

protected:
  const GridView gridView_;
  BilinearForm& bilinearForm;
  InnerProduct& innerProduct;
  OptimalTestIndexSet<BilinearForm, InnerProduct> indexSet_;
};




/** \brief The restriction of a finite element basis to a single element */
template<typename BilinearForm, typename InnerProduct>
class OptimalTestBasisLocalView
{
public:

  typedef typename BilinearForm::SolutionSpaces SolutionSpaces;
  typedef typename BilinearForm::TestSpaces EnrichedTestspaces;
  /** \brief The global FE basis that this is a view on */
  typedef OptimalTestBasis<BilinearForm,InnerProduct> GlobalBasis;
  typedef typename GlobalBasis::GridView GridView;

  /** \brief The type used for sizes */
  typedef typename GlobalBasis::size_type size_type;

  /** \brief Type used to number the degrees of freedom
   *
   * In the case of mixed finite elements this really can be a multi-index, but for a standard
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
  typedef OptimalTestBasisLeafNode<BilinearForm,InnerProduct> Tree;

  /** \brief Construct local view for a given global finite element basis */
  OptimalTestBasisLocalView(const GlobalBasis* globalBasis,
                            BilinearForm& bilinForm,
                            InnerProduct& innerProd) :
    globalBasis_(globalBasis),
    bilinearForm(bilinForm),
    innerProduct(innerProd),
    tree_(globalBasis, bilinearForm, innerProduct)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
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
   * And indeed, in the OptimalTestBasisView implementation this method does nothing.
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
    return tree_.finiteElement_->size();
  }

  /**
   * \brief Maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   *
   * The method returns k^dim, which is the number of degrees of freedom you get for cubes.
   */
  size_type maxSize() const
  {
    return boost::fusion::fold(boost::fusion::transform(bilinearForm.getSolutionSpaces(),
                                                        getLocalViewMaxSize()),
                               0, std::plus<std::size_t>());
  }

  /** \brief Return the global basis that we are a view on
   */
  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

protected:
  const GlobalBasis* globalBasis_;
  BilinearForm& bilinearForm;
  InnerProduct& innerProduct;
  const Element* element_;
  Tree tree_;
};




template<typename BilinearForm, typename InnerProduct>
class OptimalTestBasisLeafNode :
  public GridFunctionSpaceBasisLeafNodeInterface<
    typename std::tuple_element<0,typename BilinearForm::SolutionSpaces>::type::GridView::template Codim<0>::Entity,
    Dune::OptimalTestLocalFiniteElement<typename std::tuple_element<0,typename BilinearForm::SolutionSpaces>::type::GridView::ctype,
                                        double,
                                        std::tuple_element<0,typename BilinearForm::SolutionSpaces>::type::GridView::dimension,
                                        typename std::tuple_element<0,typename BilinearForm::TestSpaces>::type::LocalView::Tree::FiniteElement  >, //TODO fuer mehrere testspaces
    typename OptimalTestBasis<BilinearForm,InnerProduct>::size_type>
{
public:
  typedef typename BilinearForm::SolutionSpaces SolutionSpaces;
  typedef typename BilinearForm::TestSpaces EnrichedTestspaces;


private:
  typedef OptimalTestBasis<BilinearForm,InnerProduct> GlobalBasis;
  typedef typename std::tuple_element<0,SolutionSpaces>::type::GridView GridView;
  static const int dim = GridView::dimension;

  typedef typename GridView::template Codim<0>::Entity E;
  typedef typename std::tuple_element<0,EnrichedTestspaces>::type::LocalView::Tree::FiniteElement EnrichedFiniteElement;
  typedef typename Dune::OptimalTestLocalFiniteElement<typename GridView::ctype,double,GridView::dimension,
                                                       EnrichedFiniteElement> FE; //TODO fuer mehrere testspaces



  typedef typename GlobalBasis::size_type ST;
  typedef typename GlobalBasis::MultiIndex MI;

  typedef typename GlobalBasis::LocalView LocalView;

  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  friend LocalView;
  friend class OptimalTestLocalIndexSet<BilinearForm, InnerProduct>;

  typedef typename boost::fusion::result_of::as_vector<
             typename boost::fusion::result_of::transform<SolutionSpaces,
                                                          getLocalView>::type
             >::type SolutionLocalView;

  typedef typename boost::fusion::result_of::as_vector<
             typename boost::fusion::result_of::transform<EnrichedTestspaces,
                                                          getLocalView>::type
             >::type TestLocalView;

public:
  typedef GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST> Interface;
  typedef typename Interface::size_type size_type;
  typedef typename Interface::Element Element;
  typedef typename Interface::FiniteElement FiniteElement;

  OptimalTestBasisLeafNode(const GlobalBasis* globalBasis, BilinearForm& bilinForm, InnerProduct& innerProd):
    globalBasis_(globalBasis),
    bilinearForm(bilinForm),
    innerProduct(innerProd),
    finiteElement_(nullptr),
    enrichedTestspace_(nullptr),
    element_(nullptr),
    localViewSolution(boost::fusion::as_vector(boost::fusion::transform(bilinForm.getSolutionSpaces(),getLocalView()))),
    localViewTest(boost::fusion::as_vector(boost::fusion::transform(bilinForm.getTestSpaces(),getLocalView())))
  { }

  //! Return current element, throw if unbound
  const Element& element() const DUNE_FINAL
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const DUNE_FINAL
  {
    return *finiteElement_;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const DUNE_FINAL
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
    return finiteElement_->size();
  }

  //! Maps from subtree index set [0..subTreeSize-1] into root index set (element-local) [0..localSize-1]
  size_type localIndex(size_type i) const DUNE_FINAL
  {
    return i;
  }

  void setLocalIndex(size_type leafindex, size_type localindex) DUNE_FINAL
  { DUNE_THROW(Dune::NotImplemented, "OptimalTestBasisLeafNode does not support setLocalIndex() yet."); }

protected:

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    for_each(localViewTest, applyBind<decltype(e)>(e));
    enrichedTestspace_ = const_cast<typename std::tuple_element<0,EnrichedTestspaces>::type::LocalView::Tree::FiniteElement*>
                          (&(boost::fusion::at_c<0>(localViewTest)->tree().finiteElement()));
    for_each(localViewSolution, applyBind<decltype(e)>(e));

    //auto solutionLocalFeSize = boost::fusion::as_vector(
    //                              boost::fusion::transform(localViewSolution_, getLocalFeSize()));
    //int n = boost::fusion::fold(solutionLocalFeSize, 0, std::plus<std::size_t>());

    Matrix<FieldMatrix<double,1,1> > stiffnessMatrix;
    innerProduct.bind(localViewTest);
    innerProduct.getLocalMatrix(stiffnessMatrix);

    bilinearForm.bind(localViewTest, localViewSolution);
    bilinearForm.getLocalMatrix(coefficientMatrix);

    Cholesky<Matrix<FieldMatrix<double,1,1> > > cholesky(stiffnessMatrix);
    cholesky.apply(coefficientMatrix);

    finiteElement_ = Dune::Std::make_unique<FiniteElement>(&coefficientMatrix, enrichedTestspace_);
  }

  const GlobalBasis* globalBasis_;
  BilinearForm& bilinearForm;
  InnerProduct& innerProduct;
  std::unique_ptr<FiniteElement> finiteElement_;
  EnrichedFiniteElement* enrichedTestspace_; //TODO fuer mehrere Testspaces
  const Element* element_;
  Matrix<FieldMatrix<double,1,1> > coefficientMatrix;
  SolutionLocalView localViewSolution;
  TestLocalView localViewTest;

};

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_OPTIMALTESTBASIS_HH
