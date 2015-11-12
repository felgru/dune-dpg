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

#include <dune/localfunctions/optimaltestfunctions/optimaltest.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
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
        if (localIndexSet.size()>*index_result)
        {
          *index_found=true;
          *index_result=(localIndexSet.index(*index_result))[0];
        }
        else
        {
          *space_index+=1;
          *index_result-=localIndexSet.size();
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
    bilinearForm_(bilinForm),
    innerProduct_(innerProd),
    gridView_(std::get<0>(bilinForm.getSolutionSpaces()).gridView()),
    localViewSolution_(boost::fusion::as_vector(
                boost::fusion::transform(bilinearForm_.getSolutionSpaces(),
                                         detail::getLocalView()))),
    localViewTest_(boost::fusion::as_vector(
                boost::fusion::transform(bilinearForm_.getTestSpaces(),
                                         detail::getLocalView()))),
    geometryBufferIsSet_(false)
  {}

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



// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   OptimalTestBasisNodeFactory
//   OptimalTestBasisNodeIndexSet
//   OptimalTestBasisNode
//
// The factory allows to create the others and is the owner of possible shared
// state. These three components do _not_ depend on the global basis or index
// set and can be used without a global basis.
// *****************************************************************************

template<typename TestspaceCoefficientMatrix, std::size_t testIndex, typename ST, typename TP>
class OptimalTestBasisNode;

template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class MI, class TP, class ST>
class OptimalTestBasisNodeIndexSet;

template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class MI, class ST>
class OptimalTestBasisNodeFactory;



template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class MI, class ST>
class OptimalTestBasisNodeFactory
{
  static const int dim =
      TestspaceCoefficientMatrix::GridView::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = typename TestspaceCoefficientMatrix::GridView;
  using size_type = ST;


  template<class TP>
  using Node = OptimalTestBasisNode<TestspaceCoefficientMatrix, testIndex, size_type, TP>;

  template<class TP>
  using IndexSet = OptimalTestBasisNodeIndexSet<TestspaceCoefficientMatrix, testIndex, MI, TP, ST>;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using SizePrefix = Dune::ReservedVector<size_type, 2>;

  using TestSearchSpaces = typename TestspaceCoefficientMatrix::TestSpaces;
  using SolutionSpaces = typename TestspaceCoefficientMatrix::SolutionSpaces;


  /** \brief Constructor for a given test coefficient matrix */
  OptimalTestBasisNodeFactory(TestspaceCoefficientMatrix& testCoeffMat) :
    testspaceCoefficientMatrix_(testCoeffMat)
  {}


  void initializeIndices()
  {
    using namespace boost::fusion;
    using namespace Dune::detail;

    /* set up global offsets */
    fold(zip(globalOffsets,
             testspaceCoefficientMatrix_.bilinearForm().getSolutionSpaces()),
         (size_t)0, globalOffsetHelper());
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return testspaceCoefficientMatrix_.gridView();
  }

  template<class TP>
  Node<TP> node(const TP& tp) const
  {
    return Node<TP>{tp, testspaceCoefficientMatrix_};
  }

  template<class TP>
  IndexSet<TP> indexSet() const
  {
    return IndexSet<TP>{*this};
  }

  size_type size() const
  {
    using namespace boost::fusion;
    using namespace Dune::detail;

    return fold(transform(testspaceCoefficientMatrix_.bilinearForm()
                              .getSolutionSpaces(),
                          getSize()),
                0, std::plus<std::size_t>());
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    if (prefix.size() == 0)
      return size();
    assert(false);
  }

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return size();
  }

  size_type maxNodeSize() const
  {
    using namespace boost::fusion;
    using namespace Dune::detail;

    return fold(transform(testspaceCoefficientMatrix_.bilinearForm()
                              .getSolutionSpaces(),
                          getMaxNodeSize()),
                0, std::plus<std::size_t>());
  }

//protected:
  // TODO: store testspaceCoefficientMatrix_ by reference or by value?
  TestspaceCoefficientMatrix& testspaceCoefficientMatrix_;

  size_t globalOffsets[std::tuple_size<SolutionSpaces>::value];

};



template<typename TestspaceCoefficientMatrix, std::size_t testIndex, typename ST, typename TP>
class OptimalTestBasisNode :
  public LeafBasisNode<ST, TP>
{
public:
  using SolutionSpaces = typename TestspaceCoefficientMatrix::SolutionSpaces;
  using TestSearchSpaces = typename TestspaceCoefficientMatrix::TestSpaces;

private:
  using GV = typename TestspaceCoefficientMatrix::GridView;
  static const int dim = GV::dimension;
  // TODO: maxSize does not seem to be needed.
  // static const int maxSize = see NodeFactory::maxNodeSize();

  using Base = LeafBasisNode<ST,TP>;
  using TestSearchFiniteElement
      = typename std::tuple_element<testIndex,TestSearchSpaces>::type
                   ::LocalView::Tree::FiniteElement;

  using SolutionLocalView
        = typename boost::fusion::result_of::as_vector<
             typename boost::fusion::result_of::transform<
                         SolutionSpaces,
                         detail::getLocalView
                      >::type
             >::type;
  using TestLocalView
        = typename boost::fusion::result_of::as_vector<
             typename boost::fusion::result_of::transform<
                         TestSearchSpaces,
                         detail::getLocalView
                      >::type
             >::type;

public:

  using size_type = ST;
  using TreePath = TP;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename Dune::OptimalTestLocalFiniteElement
                                  <typename GV::ctype,
                                   double,
                                   GV::dimension,
                                   TestSearchFiniteElement>;

  OptimalTestBasisNode(const TreePath& treePath,
                       TestspaceCoefficientMatrix& testCoeffMat) :
    Base(treePath),
    testspaceCoefficientMatrix(testCoeffMat),
    finiteElement_(nullptr),
    testSearchSpace_(nullptr),
    element_(nullptr),
    localViewSolution_(boost::fusion::as_vector(
                boost::fusion::transform(testCoeffMat.bilinearForm()
                                                .getSolutionSpaces(),
                                         detail::getLocalView()))),
    localViewTest(boost::fusion::as_vector(
                boost::fusion::transform(testCoeffMat.bilinearForm()
                                                    .getTestSpaces(),
                                         detail::getLocalView())))
  {}

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const
  {
    return *finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    using namespace Dune::detail;
    using namespace boost::fusion;

    element_ = &e;
    for_each(localViewTest, applyBind<decltype(e)>(e));
    testSearchSpace_ =
        &(at_c<testIndex>(localViewTest).tree().finiteElement());
    for_each(localViewSolution_, applyBind<decltype(e)>(e));

    testspaceCoefficientMatrix.bind(e);

    size_t k = at_c<testIndex>(localViewTest).tree().finiteElement().size();
    size_t localTestSpaceOffsets[std::tuple_size<TestSearchSpaces>::value];
    fold(zip(localTestSpaceOffsets, localViewTest), (size_t)0, offsetHelper());
    size_t offset = at_c<testIndex>(localTestSpaceOffsets);

    finiteElement_ = Dune::Std::make_unique<FiniteElement>
                        (&(testspaceCoefficientMatrix.coefficientMatrix()),
                         testSearchSpace_, offset, k);
    this->setSize(finiteElement_->size());
  }

  const SolutionLocalView& localViewSolution() const
  {
    return localViewSolution_;
  }
protected:

  TestspaceCoefficientMatrix& testspaceCoefficientMatrix;
  //FiniteElementCache cache_;
  std::unique_ptr<FiniteElement> finiteElement_;
  const TestSearchFiniteElement* testSearchSpace_;
  const Element* element_;
  SolutionLocalView localViewSolution_;
  // TODO: localViewTest is only used to get the testSearchSpace_
  TestLocalView localViewTest;
};



template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class MI, class TP, class ST>
class OptimalTestBasisNodeIndexSet
{
  using GV = typename TestspaceCoefficientMatrix::GridView;
  enum {dim = GV::dimension};

public:

  using size_type = ST;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = OptimalTestBasisNodeFactory<TestspaceCoefficientMatrix, testIndex, MI, ST>;

  using Node = typename NodeFactory::template Node<TP>;

  typedef typename TestspaceCoefficientMatrix::SolutionSpaces SolutionSpaces;
  typedef typename boost::fusion::result_of::as_vector<
             typename boost::fusion::result_of::transform<
                       SolutionSpaces,
                       detail::getLocalIndexSet>::type
             >::type SolutionLocalIndexSet;

  OptimalTestBasisNodeIndexSet(const NodeFactory& nodeFactory) :
    nodeFactory_(&nodeFactory),
    solutionLocalIndexSet_(boost::fusion::as_vector(
              boost::fusion::transform(
                      nodeFactory.testspaceCoefficientMatrix_
                          .bilinearForm().getSolutionSpaces(),
                      detail::getLocalIndexSet())))
  {}

  constexpr OptimalTestBasisNodeIndexSet(const OptimalTestBasisNodeIndexSet&)
          = default;

  OptimalTestBasisNodeIndexSet(OptimalTestBasisNodeIndexSet&& indexSet)
          = default;

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Node& node)
  {
    using namespace boost::fusion;

    node_ = &node;

    for_each(zip(solutionLocalIndexSet_, node.localViewSolution()),
             make_fused_procedure(detail::bindLocalIndexSet()));
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    node_ = nullptr;
    for_each(solutionLocalIndexSet_, detail::applyUnbind());
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    return node_->finiteElement().size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  MultiIndex index(size_type i) const
  {
    using namespace Dune::detail;
    using namespace boost::fusion;

    size_t space_index=0;
    size_t index_result=i;
    bool index_found=false;

    for_each(solutionLocalIndexSet_,
             computeIndex(space_index, index_result, index_found));

    MultiIndex result;
    result[0]=(nodeFactory_->globalOffsets[space_index]+index_result);
    return result;
  }

protected:
  const NodeFactory* nodeFactory_;

  const Node* node_;

  SolutionLocalIndexSet solutionLocalIndexSet_;
};



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Nodal basis of a scalar k-th-order Lagrangean finite element space
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * \tparam TestspaceCoefficentMatrix  caches the computation of local
 *                                    optimal test spaces
 * \tparam testIndex  index of the optimal test space in the test space tuple
 */
template<typename TestspaceCoefficientMatrix, std::size_t testIndex = 0, class ST = std::size_t>
using OptimalTestBasis = DefaultGlobalBasis<OptimalTestBasisNodeFactory<TestspaceCoefficientMatrix, testIndex, FlatMultiIndex<ST>, ST> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_OPTIMALTESTBASIS_HH
