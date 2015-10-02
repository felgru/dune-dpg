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





#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/std/final.hh>

#include <dune/localfunctions/optimaltestfunctions/optimaltest.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>
#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/cholesky.hh>


namespace Dune {
namespace Functions {




//////////////////////////////////////////////////////////////////////////////
// Helper begin
//////////////////////////////////////////////////////////////////////////////


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


//////////////////////////////////////////////////////////////////////////////
// Helper end
//////////////////////////////////////////////////////////////////////////////


template <typename Geometry>
class GeometryBuffer
{
  typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

  public:

  void set(Geometry geometry)
  {
    corners=geometry.corners();
    cornerVector.resize(corners);
    for (unsigned int i=0; i<corners; i++)
    {
      cornerVector[i]=geometry.corner(i);
    }
  };

  bool isSame(Geometry geometry)
  {
    if (geometry.corners()==corners)
    {
      GlobalCoordinate difference = cornerVector[0];
      difference -= geometry.corner(0);
      double volume = geometry.volume();    //TODO evtl. schon Volumen fuer ersten Vergleich nutzen?
      for (unsigned int i=1; i<corners; i++)
      {
        for (unsigned int j=0; j<difference.size(); j++)
        {
          if (std::abs(geometry.corner(i)[j]+difference[j]-cornerVector[i][j] )/volume > 1e-10)
          {
            return false;
          }
        }
        return true;
      }
    }
    std::cout << "different number of corners" <<std::endl;
    return false;
  };

  private:
  int corners;
  std::vector<GlobalCoordinate> cornerVector;
};




template<typename BilinForm, typename InnerProd>
class TestspaceCoefficientMatrix
{

  public:
  typedef BilinForm BilinearForm;
  typedef InnerProd InnerProduct;
  typedef typename std::tuple_element<0,typename BilinearForm::SolutionSpaces>::type::GridView GridView;
  typedef typename BilinearForm::SolutionSpaces SolutionSpaces;
  typedef typename BilinearForm::TestSpaces TestSpaces;

  private:
  typedef typename boost::fusion::result_of::as_vector<
             typename boost::fusion::result_of::transform<
                         SolutionSpaces,
                         detail::getLocalView
                      >::type
             >::type SolutionLocalView;

  typedef typename boost::fusion::result_of::as_vector<
             typename boost::fusion::result_of::transform<
                         TestSpaces,
                         detail::getLocalView
                      >::type
             >::type TestLocalView;

  typedef Matrix<FieldMatrix<double,1,1> > MatrixType;

  public:
  TestspaceCoefficientMatrix(BilinearForm& bilinForm, InnerProduct& innerProd) :
    gridView_(std::get<0>(bilinForm.getSolutionSpaces()).gridView()),
    bilinearForm_(bilinForm),
    innerProduct_(innerProd),
    localViewSolution_(boost::fusion::as_vector(
                boost::fusion::transform(bilinearForm_.getSolutionSpaces(),
                                         detail::getLocalView()))),
    localViewTest_(boost::fusion::as_vector(
                boost::fusion::transform(bilinearForm_.getTestSpaces(),
                                         detail::getLocalView()))),
    geometryBufferIsSet_(false)
  {}

  /* TODO: Can probably be replaced by default move ctor with
   *       boost::fusion 1.55. */
  TestspaceCoefficientMatrix(TestspaceCoefficientMatrix&& coeffMatrix)
  : bilinearForm_(coeffMatrix.bilinearForm_),
    innerProduct_(coeffMatrix.innerProduct_),
    gridView_(coeffMatrix.gridView_),
    geometryBuffer_(std::move(coeffMatrix.geometryBuffer_)),
    geometryBufferIsSet_(coeffMatrix.geometryBufferIsSet_),
    coefficientMatrix_(std::move(coeffMatrix.coefficientMatrix_))
  {
    using namespace boost::fusion;

    copy(coeffMatrix.localViewSolution_, localViewSolution_);
    copy(coeffMatrix.localViewTest_, localViewTest_);
    for_each(coeffMatrix.localViewSolution_, detail::setToNullptr());
    for_each(coeffMatrix.localViewTest_, detail::setToNullptr());
  }

  ~TestspaceCoefficientMatrix()
  {
    for_each(localViewTest_,     detail::default_deleter());
    for_each(localViewSolution_, detail::default_deleter());
  }

  typedef decltype(std::declval<typename GridView::template Codim<0>::Entity>().geometry()) Geometry;
  typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

  void bind(const typename GridView::template Codim<0>::Entity& e)
  {
    using namespace Dune::detail;

    if (geometryBufferIsSet_ and geometryBuffer_.isSame(e.geometry()))
    {
      //std::cout <<"Old Geometry" <<std::endl;
    }
    else
    {
      geometryBuffer_.set(e.geometry());
      //std::cout <<"New Geometry" <<std::endl;
      for_each(localViewTest_, applyBind<decltype(e)>(e));
      for_each(localViewSolution_, applyBind<decltype(e)>(e));

      MatrixType stiffnessMatrix;

      innerProduct_.bind(localViewTest_);
      innerProduct_.getLocalMatrix(stiffnessMatrix);

      bilinearForm_.bind(localViewTest_, localViewSolution_);
      bilinearForm_.getLocalMatrix(coefficientMatrix_);

      Cholesky<MatrixType> cholesky(stiffnessMatrix);
      cholesky.apply(coefficientMatrix_);
      geometryBufferIsSet_ = true;
    }
  }

  MatrixType& coefficientMatrix()
  {
    return coefficientMatrix_;
  };

  BilinearForm& bilinearForm() const
  {
    return bilinearForm_;
  };

  InnerProduct& innerProduct() const
  {
    return innerProduct_;
  };

  const GridView& gridView() const
  {
    return gridView_;
  }

  private:
  BilinearForm&     bilinearForm_;
  InnerProduct&     innerProduct_;
  GridView          gridView_;
  SolutionLocalView localViewSolution_;
  TestLocalView     localViewTest_;
  GeometryBuffer<Geometry> geometryBuffer_;
  bool              geometryBufferIsSet_;
  MatrixType        coefficientMatrix_;
};


template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
class OptimalTestBasisLocalView;

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
class OptimalTestBasisLeafNode;

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
class OptimalTestIndexSet;





template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
class OptimalTestLocalIndexSet
{
  typedef typename TestspaceCoefficientMatrix::GridView GridView;
  enum {dim = GridView::dimension};

public:
  typedef typename TestspaceCoefficientMatrix::SolutionSpaces SolutionSpaces;

  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef OptimalTestBasisLocalView<TestspaceCoefficientMatrix, testIndex> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type, 1> MultiIndex;

  OptimalTestLocalIndexSet(const SolutionSpaces& solSp)
  :
  solutionBasisIndexSet_(boost::fusion::as_vector(
              boost::fusion::transform(solSp,detail::getIndexSet()))),
  solutionLocalIndexSet_(boost::fusion::as_vector(
              transform(solutionBasisIndexSet_,detail::getLocalIndexSet())))
  {}

  /* TODO: Can probably be replaced by default move ctor with
   *       boost::fusion 1.55. */
  OptimalTestLocalIndexSet(OptimalTestLocalIndexSet&& indexSet)
  : solutionBasisIndexSet_(indexSet.solutionBasisIndexSet_)
  {
    using namespace boost::fusion;

    copy(indexSet.solutionLocalIndexSet_, solutionLocalIndexSet_);
    for_each(indexSet.solutionLocalIndexSet_, detail::setToNullptr());
    localView_ = indexSet.localView_;
    indexSet.localView_ = nullptr;
  }

  ~OptimalTestLocalIndexSet()
  {
      using namespace boost::fusion;
      for_each(solutionLocalIndexSet_, detail::default_deleter());
  }

  typedef typename boost::fusion::result_of::as_vector<
             typename boost::fusion::result_of::transform<SolutionSpaces,
                                                          detail::getIndexSet
                                                         >::type
             >::type SolutionBasisIndexSet;
  typedef typename boost::fusion::result_of::as_vector<
             typename boost::fusion::result_of::transform<
                       SolutionBasisIndexSet,
                       detail::getLocalIndexSet>::type
             >::type SolutionLocalIndexSet;

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually
   * access any of its data members offers to centralize some expensive
   * setup code in the 'bind' method, which can save a lot of run-time costs.
   */
  void bind(const OptimalTestBasisLocalView<TestspaceCoefficientMatrix, testIndex>& localView)
  {
    using namespace Dune::detail;
    using namespace boost::fusion;

    localView_ = &localView;

    for_each(zip(solutionLocalIndexSet_, localView_->tree().localViewSolution),
             make_fused_procedure(bindLocalIndexSet()));
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    using namespace Dune::detail;

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
    using namespace Dune::detail;
    using namespace boost::fusion;

    /* set up global offsets*/
    size_t globalOffsets[std::tuple_size<SolutionSpaces>::value];
    fold(zip(globalOffsets, solutionBasisIndexSet_),
         (size_t)0, globalOffsetHelper());

    size_t space_index=0;
    size_t index_result=i;
    bool index_found=false;

    for_each(solutionLocalIndexSet_,
             computeIndex(space_index, index_result, index_found));

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

  const OptimalTestBasisLocalView<TestspaceCoefficientMatrix, testIndex>*
                                                                   localView_;

  SolutionBasisIndexSet solutionBasisIndexSet_;
  SolutionLocalIndexSet solutionLocalIndexSet_;
};



template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
class OptimalTestIndexSet
{
  typedef typename TestspaceCoefficientMatrix::GridView GridView;
  enum {dim = GridView::dimension};

  // Needs the mapper
  friend class OptimalTestLocalIndexSet<TestspaceCoefficientMatrix,testIndex>;

public:
  typedef typename TestspaceCoefficientMatrix::SolutionSpaces SolutionSpaces;

  typedef OptimalTestLocalIndexSet<TestspaceCoefficientMatrix, testIndex> LocalIndexSet;

  OptimalTestIndexSet(const SolutionSpaces& solSp)
  : solutionSpaces(solSp),
    gridView_(std::get<0>(solSp).gridView())
  {}

  std::size_t size() const
  {
    using namespace Dune::detail;

    return boost::fusion::fold(boost::fusion::transform(solutionSpaces,
                                                        getIndexSetSize()),
                               0, std::plus<std::size_t>());
  }

  LocalIndexSet localIndexSet() const
  {
    return LocalIndexSet(solutionSpaces);
  }

private:
  const SolutionSpaces& solutionSpaces;
  const GridView gridView_;
};



/** \brief Basis for DPG optimal test spaces
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - Grids must be 1d, 2d, or 3d
 * - 3d grids must be simplex grids
 *
 * \tparam TestspaceCoefficentMatrix  caches the computation of local
 *                                    optimal test spaces
 * \tparam testIndex  index of the optimal test space in the test space tuple
 */

template<typename TestspaceCoefficientMatrix, std::size_t testIndex = 0>
class OptimalTestBasis
: public GridViewFunctionSpaceBasis<typename TestspaceCoefficientMatrix::GridView,
                                    OptimalTestBasisLocalView<TestspaceCoefficientMatrix, testIndex>,
                                    OptimalTestIndexSet<TestspaceCoefficientMatrix, testIndex>,
                                    std::array<std::size_t, 1> >
{
public:
  /** \brief The grid view that the FE space is defined on */
  typedef typename TestspaceCoefficientMatrix::GridView GridView;

private:
  enum {dim = GridView::dimension};

public:
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef OptimalTestBasisLocalView<TestspaceCoefficientMatrix, testIndex> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type, 1> MultiIndex;   //TODO: Brauche ich das?

  typedef typename TestspaceCoefficientMatrix::SolutionSpaces SolutionSpaces;
  typedef typename TestspaceCoefficientMatrix::TestSpaces EnrichedTestspaces;

  /** \brief Constructor for a given grid view object */
  OptimalTestBasis(TestspaceCoefficientMatrix& testCoeffMat) :
    gridView_(testCoeffMat.gridView()),
    testspaceCoefficientMatrix(testCoeffMat),
    indexSet_(testCoeffMat.bilinearForm().getSolutionSpaces())
  {}

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const DUNE_FINAL
  {
    return gridView_;
  }

  OptimalTestIndexSet<TestspaceCoefficientMatrix,testIndex> indexSet() const
  {
    return indexSet_;
  }

  /** \brief Return local view for basis
   *
   */
  LocalView localView() const
  {
    return LocalView(this, testspaceCoefficientMatrix);
  }

protected:
  const GridView gridView_;
  TestspaceCoefficientMatrix& testspaceCoefficientMatrix;
  OptimalTestIndexSet<TestspaceCoefficientMatrix, testIndex> indexSet_;
};




/** \brief The restriction of a finite element basis to a single element */
template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
class OptimalTestBasisLocalView
{
public:

  typedef typename TestspaceCoefficientMatrix::SolutionSpaces SolutionSpaces;
  typedef typename TestspaceCoefficientMatrix::TestSpaces EnrichedTestspaces;
  /** \brief The global FE basis that this is a view on */
  typedef OptimalTestBasis<TestspaceCoefficientMatrix, testIndex> GlobalBasis;
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
  typedef OptimalTestBasisLeafNode<TestspaceCoefficientMatrix, testIndex> Tree;

  /** \brief Construct local view for a given global finite element basis */
  OptimalTestBasisLocalView(const GlobalBasis* globalBasis,
                            TestspaceCoefficientMatrix& testCoeffMat) :
    globalBasis_(globalBasis),
    testspaceCoefficientMatrix(testCoeffMat),
    tree_(globalBasis, testCoeffMat)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually
   * access any of its data members offers to centralize some expensive
   * setup code in the 'bind' method, which can save a lot of run-time costs.
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
    using namespace boost::fusion;
    using namespace Dune::detail;

    return fold(transform(testspaceCoefficientMatrix.bilinearForm()
                          .getSolutionSpaces(),
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
  TestspaceCoefficientMatrix& testspaceCoefficientMatrix;
  const Element* element_;
  Tree tree_;
};




template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
class OptimalTestBasisLeafNode :
  public GridFunctionSpaceBasisLeafNodeInterface<
    typename TestspaceCoefficientMatrix::GridView::template Codim<0>::Entity,
    Dune::OptimalTestLocalFiniteElement<typename TestspaceCoefficientMatrix::GridView::ctype,
                                        double,
                                        TestspaceCoefficientMatrix::GridView::dimension,
                                        typename std::tuple_element<testIndex, typename TestspaceCoefficientMatrix::BilinearForm::TestSpaces>::type::LocalView::Tree::FiniteElement  >,
    typename OptimalTestBasis<TestspaceCoefficientMatrix,testIndex>::size_type>
{
public:
  typedef typename TestspaceCoefficientMatrix::SolutionSpaces SolutionSpaces;
  typedef typename TestspaceCoefficientMatrix::TestSpaces EnrichedTestspaces;
  typedef typename TestspaceCoefficientMatrix::BilinearForm BilinearForm;
    typedef typename TestspaceCoefficientMatrix::InnerProduct InnerProduct;

private:
  typedef OptimalTestBasis<TestspaceCoefficientMatrix,testIndex> GlobalBasis;
  typedef typename TestspaceCoefficientMatrix::GridView GridView;
  enum {dim = GridView::dimension};

  typedef typename GridView::template Codim<0>::Entity E;
  typedef typename std::tuple_element<testIndex,EnrichedTestspaces>::type::LocalView::Tree::FiniteElement EnrichedFiniteElement;
  typedef typename Dune::OptimalTestLocalFiniteElement<typename GridView::ctype,double,GridView::dimension,
                                                       EnrichedFiniteElement> FE;


  typedef typename GlobalBasis::size_type ST;
  typedef typename GlobalBasis::MultiIndex MI;

  typedef typename GlobalBasis::LocalView LocalView;

//  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  friend LocalView;
  friend class OptimalTestLocalIndexSet<TestspaceCoefficientMatrix,testIndex>;

  typedef typename boost::fusion::result_of::as_vector<
             typename boost::fusion::result_of::transform<
                         SolutionSpaces,
                         detail::getLocalView
                      >::type
             >::type SolutionLocalView;
  typedef typename boost::fusion::result_of::as_vector<
             typename boost::fusion::result_of::transform<
                         EnrichedTestspaces,
                         detail::getLocalView
                      >::type
             >::type TestLocalView;

public:
  typedef GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST> Interface;
  typedef typename Interface::size_type size_type;
  typedef typename Interface::Element Element;
  typedef typename Interface::FiniteElement FiniteElement;

  OptimalTestBasisLeafNode(const GlobalBasis* globalBasis, TestspaceCoefficientMatrix& testCoeffMat):
    globalBasis_(globalBasis),
    testspaceCoefficientMatrix(testCoeffMat),
    bilinearForm(testspaceCoefficientMatrix.bilinearForm()),
    innerProduct(testspaceCoefficientMatrix.innerProduct()),
    finiteElement_(nullptr),
    enrichedTestspace_(nullptr),
    element_(nullptr),
    localViewSolution(boost::fusion::as_vector(
                boost::fusion::transform(bilinearForm.getSolutionSpaces(),
                                         detail::getLocalView()))),
    localViewTest(boost::fusion::as_vector(
                boost::fusion::transform(bilinearForm.getTestSpaces(),
                                         detail::getLocalView())))
  { }

  /* TODO: Can probably be replaced by default move ctor with
   *       boost::fusion 1.55. */
  OptimalTestBasisLeafNode(OptimalTestBasisLeafNode&& ln)
  : globalBasis_(std::move(ln.globalBasis_)),
    testspaceCoefficientMatrix(ln.testspaceCoefficientMatrix),
    bilinearForm(ln.bilinearForm),
    innerProduct(ln.innerProduct),
    finiteElement_(std::move(ln.finiteElement_)),
    enrichedTestspace_(std::move(ln.enrichedTestspace_)),
    element_(std::move(ln.element_))
  {
    using namespace boost::fusion;

    copy(ln.localViewSolution, localViewSolution);
    copy(ln.localViewTest, localViewTest);
    for_each(ln.localViewSolution, detail::setToNullptr());
    for_each(ln.localViewTest, detail::setToNullptr());
  }

  ~OptimalTestBasisLeafNode()
  {
    for_each(localViewTest,     detail::default_deleter());
    for_each(localViewSolution, detail::default_deleter());
  }

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
    using namespace Dune::detail;
    using namespace boost::fusion;

    element_ = &e;
    for_each(localViewTest, applyBind<decltype(e)>(e));
    enrichedTestspace_ =
        const_cast<typename std::tuple_element<testIndex,EnrichedTestspaces>
                   ::type::LocalView::Tree::FiniteElement*>
                  (&(at_c<testIndex>(localViewTest)->tree().finiteElement()));
    for_each(localViewSolution, applyBind<decltype(e)>(e));

    testspaceCoefficientMatrix.bind(e);

    // coefficientMatrix = testspaceCoefficientMatrix.coefficientMatrix();

    size_t k = at_c<testIndex>(localViewTest)->tree().finiteElement().size();
    size_t localTestSpaceOffsets[std::tuple_size<EnrichedTestspaces>::value];
    fold(zip(localTestSpaceOffsets, localViewTest), (size_t)0, offsetHelper());
    size_t offset = at_c<testIndex>(localTestSpaceOffsets);

    finiteElement_ = Dune::Std::make_unique<FiniteElement>
                        (&(testspaceCoefficientMatrix.coefficientMatrix()),
                         enrichedTestspace_, offset, k);

  }

  const GlobalBasis* globalBasis_;
  TestspaceCoefficientMatrix& testspaceCoefficientMatrix;
  BilinearForm& bilinearForm;
  InnerProduct& innerProduct;
  std::unique_ptr<FiniteElement> finiteElement_;
  EnrichedFiniteElement* enrichedTestspace_;
  const Element* element_;
//  Matrix<FieldMatrix<double,1,1> > coefficientMatrix;
  SolutionLocalView localViewSolution;
  TestLocalView localViewTest;
};

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_OPTIMALTESTBASIS_HH
