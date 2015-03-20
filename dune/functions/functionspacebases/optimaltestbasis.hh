// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_OPTIMALTESTBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_OPTIMALTESTBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
#include <dune/common/std/final.hh>
#else
 #ifndef DUNE_FINAL
  #define DUNE_FINAL
 #endif
#endif

//#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/localfunctions/optimaltestfunctions/optimaltest.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>


namespace Dune {
namespace Functions {

template<typename GV, typename EnrichedTestspace, typename SolutionSpace>
class OptimalTestBasisLocalView;

template<typename GV, typename EnrichedTestspace, typename SolutionSpace>
class OptimalTestBasisLeafNode;

template<typename GV, typename EnrichedTestspace, typename SolutionSpace>
class OptimalTestIndexSet;

template<typename GV, typename EnrichedTestspace, typename SolutionSpace>
class OptimalTestLocalIndexSet
{
  enum {dim = GV::dimension};

public:
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef OptimalTestBasisLocalView<GV, EnrichedTestspace, SolutionSpace> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type, 1> MultiIndex;

  OptimalTestLocalIndexSet(const OptimalTestIndexSet<GV, EnrichedTestspace, SolutionSpace> & indexSet)
  : basisIndexSet_(indexSet),
  solutionLocalIndexSet_(indexSet.solutionSpace.indexSet())
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const OptimalTestBasisLocalView<GV,EnrichedTestspace, SolutionSpace>& localView)
  {
    localView_ = &localView;
    solutionLocalIndexSet_.bind(localView_->tree().localViewSolution);
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    localView_ = nullptr;
    solutionLocalIndexSet_.unbind();
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
    return localView_->tree().finiteElement_->size();
#else
    return localView_->tree().finiteElement_->localBasis().size();
#endif
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  const MultiIndex index(size_type i) const
  {
    return solutionLocalIndexSet_.index(i);
    DUNE_THROW(Dune::NotImplemented, "Grid contains elements not supported for the PQKNodalBasis");
  }

  /** \brief Return the local view that we are attached to
   */
  const LocalView& localView() const
  {
    return *localView_;
  }

  const OptimalTestBasisLocalView<GV,EnrichedTestspace, SolutionSpace>* localView_;

  const OptimalTestIndexSet<GV,EnrichedTestspace, SolutionSpace> basisIndexSet_;

  typename OptimalTestIndexSet<GV, EnrichedTestspace, SolutionSpace>::SolutionIndexSet::LocalIndexSet solutionLocalIndexSet_;
};





template<typename GV, typename EnrichedTestspace, typename SolutionSpace>
class OptimalTestIndexSet
{
  static const int dim = GV::dimension;

  // Needs the mapper
  friend class OptimalTestLocalIndexSet<GV, EnrichedTestspace, SolutionSpace>;

public:

  typedef OptimalTestLocalIndexSet<GV, EnrichedTestspace, SolutionSpace> LocalIndexSet;
  typedef decltype(std::declval<SolutionSpace>().indexSet()) SolutionIndexSet;

  OptimalTestIndexSet(const GV& gridView)
  : gridView_(gridView),
  solutionSpace(gridView)
  {}

  std::size_t size() const
  {
    return solutionSpace.indexSet().size();
    //DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
  }

  LocalIndexSet localIndexSet() const
  {
    return LocalIndexSet(*this);
  }

private:
  const GV gridView_;
  SolutionSpace solutionSpace;
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
template<typename GV, typename EnrichedTestspace, typename SolutionSpace>
class OptimalTestBasis
: public GridViewFunctionSpaceBasis<GV,
                                    OptimalTestBasisLocalView<GV,EnrichedTestspace,SolutionSpace>,
                                    OptimalTestIndexSet<GV,EnrichedTestspace,SolutionSpace>,
                                    std::array<std::size_t, 1> >
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  typedef GV GridView;
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef OptimalTestBasisLocalView<GV,EnrichedTestspace,SolutionSpace> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  //typedef std::array<size_type, 1> MultiIndex;   //TODO: Brauche ich das?

  /** \brief Constructor for a given grid view object */
  OptimalTestBasis(const GridView& gv) :
    gridView_(gv),
    indexSet_(gv)
  {}

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const DUNE_FINAL
  {
    return gridView_;
  }

  OptimalTestIndexSet<GV, EnrichedTestspace, SolutionSpace> indexSet() const
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

protected:
  const GridView gridView_;

  OptimalTestIndexSet<GV, EnrichedTestspace, SolutionSpace> indexSet_;
};


/** \brief The restriction of a finite element basis to a single element */
template<typename GV, typename EnrichedTestspace, typename SolutionSpace>
class OptimalTestBasisLocalView
{
public:
  /** \brief The global FE basis that this is a view on */
  typedef OptimalTestBasis<GV,EnrichedTestspace,SolutionSpace> GlobalBasis;
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
  typedef OptimalTestBasisLeafNode<GV,EnrichedTestspace,SolutionSpace> Tree;

  /** \brief Construct local view for a given global finite element basis */
  OptimalTestBasisLocalView(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    tree_(globalBasis)
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
   * The method returns k^dim, which is the number of degrees of freedom you get for cubes.
   */
  size_type maxSize() const
  {
    return 4; //StaticPower<(k+1),GV::dimension>::power; //TODO TrialSpace maxSize()
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


template<typename GV, typename EnrichedTestspace, typename SolutionSpace>
class OptimalTestBasisLeafNode :
  public GridFunctionSpaceBasisLeafNodeInterface<
    typename GV::template Codim<0>::Entity,
    Dune::OptimalTestLocalFiniteElement<typename GV::ctype,double,GV::dimension, EnrichedTestspace>,
    typename OptimalTestBasis<GV,EnrichedTestspace,SolutionSpace>::size_type>
{
  typedef OptimalTestBasis<GV,EnrichedTestspace,SolutionSpace> GlobalBasis;
  static const int dim = GV::dimension;

  typedef typename GV::template Codim<0>::Entity E;
  typedef typename Dune::OptimalTestLocalFiniteElement<typename GV::ctype,double,GV::dimension, EnrichedTestspace> FE;
  typedef typename GlobalBasis::size_type ST;
//  typedef typename GlobalBasis::MultiIndex MI;

  typedef typename GlobalBasis::LocalView LocalView;

  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  friend LocalView;
  friend class OptimalTestLocalIndexSet<GV, EnrichedTestspace, SolutionSpace>;

public:
  typedef GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST> Interface;
  typedef typename Interface::size_type size_type;
  typedef typename Interface::Element Element;
  typedef typename Interface::FiniteElement FiniteElement;

  OptimalTestBasisLeafNode(const GlobalBasis* globalBasis):
    globalBasis_(globalBasis),
    finiteElement_(nullptr),
    element_(nullptr),
    solutionSpace(globalBasis->gridView()),
    localViewSolution(&solutionSpace)
  {
    coefficientMatrix.setBuildMode(MatrixType::implicit);
    coefficientMatrix.setImplicitBuildModeParameters(enrichedTestspace.size(),0);
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

protected:

  //! Bind to element.
  void bind(const Element& e)
  {
    if (finiteElement_ != nullptr)
    {
      delete finiteElement_;
    }
    element_ = &e;
    localViewSolution.bind(e);
    int n = localViewSolution.tree().finiteElement().size();
    int m=enrichedTestspace.size();
    coefficientMatrix.setSize(n,m);

    computeCoefficientMatrix();

    /*for(unsigned int i=0; i<n; i++)
    {
      coefficientMatrix.entry(i,i)=1;
    }
    coefficientMatrix.entry(n-1,m-1)=1;
    coefficientMatrix.entry(0,2)=1;*/
    localKeyList_.resize(n);
    printmatrix(std::cout , coefficientMatrix, "coefficientMatrix", "--");
    for (unsigned int i=0; i<n; i++)
    {
      localKeyList_[i] = localViewSolution.tree().finiteElement().localCoefficients().localKey(i);  //TODO: Funktionier nur, wenn nur eine lokale Basis
    }
    finiteElement_ = new FiniteElement(&coefficientMatrix, &localKeyList_); //TODO ordentlich wieder freigeben
  }

  const GlobalBasis* globalBasis_;
  FiniteElement* finiteElement_;
  SolutionSpace solutionSpace;
  EnrichedTestspace enrichedTestspace;
  std::vector<LocalKey> localKeyList_;
  const Element* element_;
  MatrixType coefficientMatrix;
  typename SolutionSpace::LocalView localViewSolution;

private:
  void computeCoefficientMatrix()
  {
    FieldVector<double, 2> beta={1,1};        //TODO: Constants should be set somewhere else!
    double c=0;

    /////////////////////////////////////////////////////////
    //   Get local finite element for solution space,
    //   local finite element for enriched test space is enrichetTestspace
    /////////////////////////////////////////////////////////

    decltype(localViewSolution.tree().finiteElement()) solutionLocalFiniteElement=localViewSolution.tree().finiteElement();
    int n = solutionLocalFiniteElement.size();
    int m = enrichedTestspace.size();

    /////////////////////////////////////////////////////////
    //   Stiffness matrix and right hand side vector
    /////////////////////////////////////////////////////////

    typedef BlockVector<FieldVector<double,1> > VectorType;

    std::vector<VectorType> rhs;
    MatrixType stiffnessMatrix;

    stiffnessMatrix.setBuildMode(MatrixType::implicit);
    stiffnessMatrix.setImplicitBuildModeParameters(m,0);
    stiffnessMatrix.setSize(m,m);

    rhs.resize(n);
    for (unsigned int i=0; i<n; ++i)
    {
      rhs[i].resize(m);
      rhs[i] = 0;
    }

    /////////////////////////////////////////////////////////
    //   Get geometry
    /////////////////////////////////////////////////////////

    const int dim = Element::dimension;
    auto geometry = element_->geometry();

    for (unsigned int i=0; i<n; ++i)
    {
      rhs[i].resize(m);
      rhs[i] = 0;
    }

    /////////////////////////////////////////////////////////
    //   Get a quadrature rule
    /////////////////////////////////////////////////////////

    int order = 2*(dim*enrichedTestspace.localBasis().order()-1);
    const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element_->type(), order);

    /////////////////////////////////////////////////////////
    //   Loop over all quadrature points
    /////////////////////////////////////////////////////////
    for (size_type pt=0; pt < quad.size(); pt++)
    {
      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();

      // The transposed inverse Jacobian of the map from the reference element to the element
      const auto& jacobian = geometry.jacobianInverseTransposed(quadPos);

      // The multiplicative factor in the integral transformation formula
      const double integrationElement = geometry.integrationElement(quadPos);

      ////////////////////////////
      // Test Functions:
      ////////////////////////////
      // values of the shape functions
      std::vector<FieldVector<double,1> > valuesTest;
      enrichedTestspace.localBasis().evaluateFunction(quadPos, valuesTest);

      // The gradients of the shape functions on the reference element
      std::vector<FieldMatrix<double,1,dim> > referenceGradientsTest;
      enrichedTestspace.localBasis().evaluateJacobian(quadPos, referenceGradientsTest);

      // Compute the shape function gradients on the real element
      std::vector<FieldVector<double,dim> > gradientsTest(referenceGradientsTest.size());
      for (size_t i=0; i<gradientsTest.size(); i++)
        jacobian.mv(referenceGradientsTest[i][0], gradientsTest[i]);
      ////////////////////////////
      // Interior Trial Functions
      ////////////////////////////
      // values of the shape functions
      std::vector<FieldVector<double,1> > valuesInterior;
      solutionLocalFiniteElement.localBasis().evaluateFunction(quadPos, valuesInterior);

      // Compute the actual matrix entries
      for (size_t i=0; i<m; i++)
      {
        for (size_t j=0; j<m; j++ )
        {
          stiffnessMatrix.entry(i,j) += (valuesTest[i] * valuesTest[j]) * quad[pt].weight() * integrationElement;
          stiffnessMatrix.entry(i,j) += (beta*gradientsTest[i]) * (beta*gradientsTest[j]) * quad[pt].weight() * integrationElement;
        }
        for (size_t j=0; j<n; j++)
        {
          rhs[j][i] += (valuesTest[i] * valuesInterior[j])* c * quad[pt].weight() * integrationElement;
          rhs[j][i] += (beta*gradientsTest[i])*(-1.0) * valuesInterior[j] * quad[pt].weight() * integrationElement;
        }
      }
    }

    UMFPack<MatrixType> umfPack(stiffnessMatrix, 0);
    InverseOperatorResult statistics;

    for (unsigned int i = 0; i<n; i++)
    {
      VectorType solution(m);
      solution = 0;

      ////////////////////////////
      //   Compute solution
      ////////////////////////////

      //UMFPack<MatrixType> umfPack(stiffnessMatrix, 2);      \\TODO: gucken, ob das so geht
      //InverseOperatorResult statistics;
      umfPack.apply(solution, rhs[i], statistics);

      for (unsigned int j = 0; j<m; j++)
      {
        coefficientMatrix.entry(i,j)=solution[j];
      }
    }
/*
    // Technicality:  turn the matrix into a linear operator
    MatrixAdapter<MatrixType,VectorType,VectorType> op(stiffnessMatrix);

    // Sequential incomplete LU decomposition as the preconditioner
    SeqILU0<MatrixType,VectorType,VectorType> ilu0(stiffnessMatrix,1.0);

    // Preconditioned conjugate-gradient solver
    CGSolver<VectorType> cg(op,
                            ilu0, // preconditioner
                            1e-4, // desired residual reduction factor
                            50,   // maximum number of iterations
                            2);   // verbosity of the solver

    // Object storing some statistics about the solving process
    InverseOperatorResult statistics;

    // Solve!
    cg.apply(x, rhs, statistics);*/
  }

};

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_OPTIMALTESTBASIS_HH
