using namespace Dune;

// HACK: THE FOLLOWING CLASS SHOULD REALLY BE AVAILABLE FROM THE DUNE-FUNCTIONS MODULE
//  However, currently it isn't, so I add it here to get a preliminary working program.


/** \brief A VTK basis grid function.
 *
 *  This function "evaluates" by evaluating the global basis and interpolating the function values.
 *  \tparam Basis   The global basis.
 *  \tparam CoefficientType The vector type for the coefficients.
*/
template<class Basis, class CoefficientType>
class VTKBasisGridFunction : public Dune::VTKFunction<typename Basis::GridView>
{
  typedef Dune::VTKFunction<typename Basis::GridView> Base;
  typedef typename  CoefficientType::value_type RangeType;
public:
  typedef typename Base::Entity Entity;
  typedef typename Base::ctype ctype;
  using Base::dim;

  /** \brief Construct from given global basis, coefficient vector and name.
   *
   *  \param basis    The global basis.
   *  \param v    A corresponding vector of coefficients.
   *  \param s    A name of the function.
   */
  VTKBasisGridFunction(const Basis &basis, const CoefficientType &v, const std::string &s) :
      basis_(basis),
      coeffs_(v),
      s_( s )
  {
    if (v.size() !=basis_.subIndexCount())
      DUNE_THROW(Dune::IOError, "VTKGridFunction: Coefficient vector is not matching the basis");
  }

  /** \brief Get the number of components the function has. */
  virtual int ncomps () const
  {
    return CoefficientType::value_type::dimension;
  }

  /** \brief Locally evaluate a component of the function.
   *
   *  \param comp The component to evaluate.
   *  \param e    The element the local coordinates are taken from.
   *  \param xi   The local coordinates where to evaluate the function.
   */
  virtual double evaluate (int comp, const Entity &e,
                           const Dune::FieldVector<ctype,dim> &xi) const
  {
    typename Basis::LocalView localView(&basis_);
    localView.bind(e);

    std::vector<RangeType> shapeFunctionValues;
    auto& basis = localView.tree().finiteElement().localBasis();
    basis.evaluateFunction(xi,shapeFunctionValues);
    RangeType r = 0;
    for (size_t i = 0; i < basis.size(); ++i)
      r += coeffs_[localView.tree().globalIndex(i)[0]] * shapeFunctionValues[i];

    return r[comp];
  }

  /** \brief Get the name of that function. */
  virtual std::string name () const
  {
    return s_;
  }

  /** \brief Destructor. */
  virtual ~VTKBasisGridFunction() {}

private:
  const Basis &basis_;
  const CoefficientType &coeffs_;
  std::string s_;
};


/**
 *  \brief Virtual interface class for grid functions
 *
 *  A general grid function that might be defined only on parts of the underlying grid (isDefinedOn).
 *  If you have a grid function defined on a GridView derive from VirtualGridViewFunction instead.
 *
 *  \tparam GridType the underlying GridType
 *  \tparam RType the type of the range space
 *
 */
template <class GridType, class RType>
class VirtualGridFunction :
  public Dune::VirtualFunction<typename GridType::template Codim<0>::Geometry::GlobalCoordinate, RType>
{
  protected:

    typedef Dune::VirtualFunction<typename GridType::template Codim<0>::Geometry::GlobalCoordinate, RType> BaseType;
    typedef VirtualGridFunction<GridType, RType> ThisType;

  public:

    typedef typename GridType::template Codim<0>::Geometry::LocalCoordinate LocalDomainType;
    typedef typename BaseType::DomainType DomainType;
    typedef typename BaseType::RangeType RangeType;


    //! The grid type
    typedef GridType Grid;

    //! Element type of underlying grid
    typedef typename Grid::template Codim<0>::Entity Element;

    /** \brief Constructor
     *
     *  \param grid the underlying grid
     */
    VirtualGridFunction(const Grid& grid):
      grid_(&grid)
    {}

    /**
     * \brief Function evaluation in local coordinates.
     *
     * \param e Evaluate in local coordinates with respect to this element.
     * \param x Argument for function evaluation in local coordinates.
     * \param y Result of function evaluation.
     */
    virtual void evaluateLocal(const Element& e, const LocalDomainType& x, RangeType& y) const = 0;


    /**
     * \brief Check whether local evaluation is possible
     *
     * \param e Return true if local evaluation is possible on this element.
     */
    virtual bool isDefinedOn(const Element& e) const = 0;

    virtual ~VirtualGridFunction() {}

protected:
    //! the underlying grid
    const Grid* grid_;

}; // end of VirtualGridFunction class





/**
 * \brief Virtual interface class for grid functions associated to grid views
 *
 */
template <class GridViewType, class RType>
class VirtualGridViewFunction :
  public VirtualGridFunction<typename GridViewType::Grid, RType>
{
  protected:
    typedef VirtualGridFunction<typename GridViewType::Grid, RType> BaseType;

  public:

    typedef typename BaseType::LocalDomainType LocalDomainType;
    typedef typename BaseType::DomainType DomainType;
    typedef typename BaseType::RangeType RangeType;
    typedef typename BaseType::DerivativeType DerivativeType;

    typedef GridViewType GridView;
    typedef typename BaseType::Grid Grid;

    //! Element type of underlying grid
    typedef typename Grid::template Codim<0>::Entity Element;

    /**
     *  \brief Constructor
     *
     *  We need to know the underlying grid view
     */
    VirtualGridViewFunction(const GridView& gridView) :
      BaseType(gridView.grid()),
      gridView_(gridView)
    {}

    /**
     * \brief Function evaluation in local coordinates.
     *
     * \param e Evaluate in local coordinates with respect to this element.
     * \param x Argument for function evaluation in local coordinates.
     * \param y Result of function evaluation.
     */
    virtual void evaluateLocal(const Element& e, const LocalDomainType& x, RangeType& y) const = 0;


    virtual ~VirtualGridViewFunction() {}

    /**
     * \brief Export underlying grid view
     */
    const GridView& gridView() const
    {
      return gridView_;
    }


  protected:

    const GridView gridView_;

}; // end of VirtualGridViewFunction class


/** \brief Cached evaluation of local basis functions
 *
 * This class allows to evaluate single local basis functions efficiently
 * by providing a cache for the evaluation of all functions.
 *
 * \tparam Imp Implementation providing the evaluateAll method
 * \tparam FunctionBaseClass The base class for this wrapper
 */
template <class Imp, class AllRangeType, class FunctionBaseClass>
class CachedComponentWrapper :
    public FunctionBaseClass
{
    public:
        typedef typename FunctionBaseClass::DomainType DomainType;
        typedef typename FunctionBaseClass::RangeType RangeType;

    private:
        typedef DomainType Position;
        struct DomainCmp
        {
            bool operator() (const DomainType& a, const DomainType& b) const
            {
                for(int i=0; i<DomainType::dimension; ++i)
                {
                    if (a[i]<b[i])
                        return true;
                    if (a[i]>b[i])
                        return false;
                }
                return false;
            }
        };
        typedef typename std::map<DomainType, AllRangeType, DomainCmp> EvaluationCache;

    public:

        CachedComponentWrapper(const int i=0) :
            i_(i)
        {}

        /** \brief Select local basis function to evaluate.
         *
         * \param i Index of the local basis function within the local basis to evaluate next.
         */
        void setIndex(int i)
        {
            i_ = i;
        }

        /** \brief Evaluate local basis function selected last using setIndex
         *
         * \param x Local coordinates with repect to entity given in constructor.
         * \param y The result is stored here.
         */
        void evaluate(const DomainType& x, RangeType& y) const
        {
            typename EvaluationCache::const_iterator it = evalCache_.find(x);
            if (it != evalCache_.end())
                y = it->second[i_];
            else
            {
                AllRangeType& yy = evalCache_[x];
                asImp().evaluateAll(x, yy);
                y = yy[i_];
            }
        }

    protected:

        const Imp& asImp () const
        {
            return static_cast<const Imp &> (*this);
        }

        int i_;

        mutable EvaluationCache evalCache_;
};


template<class F, class Grid, class Base>
class LocalFunctionComponentWrapper :
    public CachedComponentWrapper<
        LocalFunctionComponentWrapper<F, Grid, Base>,
        typename F::RangeType,
        Base>
{
    private:
        CachedComponentWrapper<LocalFunctionComponentWrapper<F, Grid, Base>, typename F::RangeType, Base> BaseClass;
        typedef typename Grid::template Codim<0>::Entity Element;
        typedef typename Element::Geometry Geometry;
    public:
        typedef typename F::RangeType AllRangeType;

        typedef typename Base::DomainType DomainType;
        typedef typename Base::RangeType RangeType;

        LocalFunctionComponentWrapper(const F& f, const Element& e, int comp) :
            BaseClass(comp),
            f_(f),
            e_(e),
            geometry_(e_.geometry())
        {}

        void evaluateAll(const DomainType& x, AllRangeType& y) const
        {
            typedef VirtualGridFunction<Grid, AllRangeType> GridFunction;

            const GridFunction* gf = dynamic_cast<const GridFunction*>(&f_);
            if (gf and gf->isDefinedOn(e_))
                gf->evaluateLocal(e_, x, y);
            else
                f_.evaluate(geometry_.global(x), y);
        }


    protected:
        const F& f_;
        const Element& e_;
        const Geometry geometry_;
};


/**
 * \brief Interpolate given function in discrete function space
 *
 * \param basis Global function space basis of discrete function space
 * \param coeff Coefficient vector to represent the interpolation
 * \param f Function to interpolate
 */
template <class B, class C, class F>
void interpolate(const B& basis, C& coeff, const F& f)
{
  typedef typename B::GridView GridView;
  typedef typename GridView::Grid Grid;

  typedef typename B::LocalView::Tree::FiniteElement LocalFiniteElement;

  typedef typename Dune::LocalFiniteElementFunctionBase<LocalFiniteElement>::type FunctionBaseClass;
  typedef LocalFunctionComponentWrapper<F, Grid, FunctionBaseClass> LocalWrapper;

  const GridView& gridview = basis.gridView();

  coeff.resize(basis.subIndexCount());

  typename Dune::BitSetVector<1> processed(basis.subIndexCount(), false);

  std::vector<typename LocalWrapper::RangeType> interpolationValues;

  typename B::LocalView localView(&basis);

  auto it = gridview.template begin<0>();
  auto end = gridview.template end<0>();
  for(; it != end; ++it)
  {
    localView.bind(*it);

    const auto& fe = localView.tree().finiteElement();

    // check if all components have already been processed
    bool allProcessed = true;
    for (size_t i=0; i<fe.localBasis().size(); ++i)
      allProcessed = allProcessed and processed[localView.tree().globalIndex(i)[0]][0];

    if (not(allProcessed))
    {
      LocalWrapper fjLocal(f, *it, 0);

      //size_t lastComponent = Components<typename F::RangeType>::size(coeff[0])-1;
      size_t lastComponent = coeff[0].size()-1;
      for (size_t j=0; j<=lastComponent; ++j)
      {
        fjLocal.setIndex(j);
        fe.localInterpolation().interpolate(fjLocal, interpolationValues);

        for (size_t i=0; i<fe.localBasis().size(); ++i)
        {
          size_t index = localView.tree().globalIndex(i)[0];

          // check if DOF is unprocessed
          if (not(processed[index][0]))
          {
            coeff[index][j] = interpolationValues[i];
            if (j==lastComponent)
              processed[index][0] = true;
          }
        }
      }
    }
  }
}
