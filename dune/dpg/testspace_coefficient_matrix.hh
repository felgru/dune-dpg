// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_TESTSPACECOEFFICIENTMATRIX_HH
#define DUNE_TESTSPACECOEFFICIENTMATRIX_HH

#include <cassert>
#include <memory>
#include <tuple>
#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dune/dpg/assemble_helper.hh>

#if DUNE_DPG_USE_LEAST_SQUARES_INSTEAD_OF_CHOLESKY
#  include <dune/dpg/leastsquares.hh>
#else
#  include <dune/dpg/cholesky.hh>
#endif

namespace Dune {

template <typename Geometry>
class GeometryComparable
{
  using GlobalCoordinate = typename Geometry::GlobalCoordinate;

  public:
  void set(const Geometry& geometry)
  {
    corners=geometry.corners();
    cornerVector.resize(corners);
    for (unsigned int i=0; i<corners; i++)
    {
      cornerVector[i]=geometry.corner(i);
      // First, I thought it would make sense to order the corners in some way
      // but then I realized, that two triangles where one is just a shift of
      // the other result in different coefficients of the optimal testspace
      // if the corners are ordered in a different way because the
      // basis functions on the physical domain are different
#if 0
      unsigned int j=0;
      while (j<i and (cornerVector[j][0]<geometry.corner(i)[0]
                      or (cornerVector[j][0]==geometry.corner(i)[0]
                          and cornerVector[j][1]<=geometry.corner(i)[1])))
      {
        j++;
      }
      for (unsigned int k=i; k>j; k--)
      {
        cornerVector[k]=cornerVector[k-1];
      }
      cornerVector[j]=geometry.corner(i);
#endif
    }
  }

  unsigned int getNumCorners() const
  {
    return corners;
  }

  GlobalCoordinate getCorner(unsigned int i) const
  {
    assert(i<corners);
    return cornerVector[i];
  }

  bool operator< (const GeometryComparable<Geometry>& geometryComparable) const
  {
    if (corners < geometryComparable.getNumCorners())
    {
      return true;
    }
    if (corners > geometryComparable.getNumCorners())
    {
      return false;
    }
    for (unsigned int j=0; j<2; j++)    //TODO works only for 2d
    {
      for (unsigned int i=1; i<corners; i++)
      {
        if ((cornerVector[i][j]-cornerVector[0][j])
             <(geometryComparable.getCorner(i)[j]
               -geometryComparable.getCorner(0)[j]
               -((std::abs(geometryComparable.getCorner(i)[j]
                          -geometryComparable.getCorner(0)[j])
                  +std::abs(cornerVector[i][j]
                           -cornerVector[0][j])
                 )*10e-10)
              ))
        {
          return true;
        }
        if ((cornerVector[i][j]-cornerVector[0][j])
             >(geometryComparable.getCorner(i)[j]
               -geometryComparable.getCorner(0)[j]
               +((std::abs(geometryComparable.getCorner(i)[j]
                          -geometryComparable.getCorner(0)[j])
                  +std::abs(cornerVector[i][j]
                           -cornerVector[0][j])
                )*10e-10)
              ))
        {
          return false;
        }
      }
    }
    return false;
  }

  private:
  unsigned int                  corners;
  std::vector<GlobalCoordinate> cornerVector;
};


template <typename MatrixType>
class CoefficientMatrices
{
public:
  const MatrixType& coefficientMatrix() const
  {
    return coefficientMatrix_;
  }

  const MatrixType& systemMatrix() const
  {
    return systemMatrix_;
  }

  template<typename Entity, typename BilinearForm, typename InnerProduct,
           typename SolutionLocalViews, typename TestLocalViews>
  void set(const Entity&        e,
           BilinearForm&        bilinearForm,
           InnerProduct&        innerProduct,
           SolutionLocalViews&  localViewsSolution,
           TestLocalViews&      localViewsTest)
  {
    using namespace Dune::detail;

    bindLocalViews(localViewsTest, e);
    bindLocalViews(localViewsSolution, e);

    MatrixType stiffnessMatrix;

    innerProduct.bind(localViewsTest);
    innerProduct.getLocalMatrix(stiffnessMatrix);

    MatrixType bilinearMatrix;

    bilinearForm.bind(localViewsTest, localViewsSolution);
    bilinearForm.getLocalMatrix(bilinearMatrix);

#if DUNE_DPG_USE_LEAST_SQUARES_INSTEAD_OF_CHOLESKY
    // Solve a least-squares problem for coefficientMatrix_
    // to make sure that it is well defined even when the Gramian of
    // the inner product is severely ill conditioned.
    solveLeastSquares(stiffnessMatrix, bilinearMatrix, coefficientMatrix_);
#else
    {
      Cholesky<MatrixType> cholesky(stiffnessMatrix);
      cholesky.apply(bilinearMatrix, coefficientMatrix_);
    }
#endif

    const unsigned int m = coefficientMatrix_.M();
    systemMatrix_.setSize(m, m);

    for (unsigned int i=0; i<m; i++)
    {
      for (unsigned int j=0; j<m; j++)
      {
        systemMatrix_[i][j]=0;
        for (unsigned int k=0; k<coefficientMatrix_.N(); k++)
        {
          systemMatrix_[i][j]+=(bilinearMatrix[k][i]*coefficientMatrix_[k][j]);
        }
      }
    }
  }

private:
  MatrixType            coefficientMatrix_;     // G^{-1}B
  MatrixType            systemMatrix_;          // B^TG^{-1}B
};


template<typename BilinForm, typename InnerProd>
class UnbufferedTestspaceCoefficientMatrix
{
public:
  using BilinearForm = BilinForm;
  using InnerProduct = InnerProd;
  using SolutionSpaces = typename BilinearForm::SolutionSpaces;
  using TestSpaces     = typename BilinearForm::TestSpaces;
  using GridView = typename std::tuple_element_t<0,SolutionSpaces>::GridView;
  using Entity = typename GridView::template Codim<0>::Entity;

private:
  using SolutionLocalViews = detail::getLocalViews_t<SolutionSpaces>;
  using TestLocalViews     = detail::getLocalViews_t<TestSpaces>;

  using MatrixType = Matrix<FieldMatrix<double,1,1>>;

public:
  UnbufferedTestspaceCoefficientMatrix(BilinearForm& bilinForm,
                                       InnerProduct& innerProd) :
    bilinearForm_(bilinForm),
    innerProduct_(innerProd),
    localViewsSolution_(detail::getLocalViews(
                          *bilinearForm_.getSolutionSpaces())),
    localViewsTest_(detail::getLocalViews(*bilinearForm_.getTestSpaces()))
  {}

  void bind(const Entity& e)
  {
    coefMatrices_.set(e,
                      bilinearForm_,
                      innerProduct_,
                      localViewsSolution_,
                      localViewsTest_);
  }

  const MatrixType& coefficientMatrix() const
  {
    return coefMatrices_.coefficientMatrix();
  }

  const MatrixType& systemMatrix() const
  {
    return coefMatrices_.systemMatrix();
  }

  const BilinearForm& bilinearForm() const
  {
    return bilinearForm_;
  }

  const InnerProduct& innerProduct() const
  {
    return innerProduct_;
  }

  const GridView& gridView() const
  {
    return std::get<0>(*bilinearForm_.getSolutionSpaces()).gridView();
  }

  private:
  BilinearForm&            bilinearForm_;
  InnerProduct&            innerProduct_;
  SolutionLocalViews       localViewsSolution_;
  TestLocalViews           localViewsTest_;
  CoefficientMatrices<MatrixType> coefMatrices_;
};



template<typename Geometry>
class GeometryBuffer
{
  using MatrixType = Matrix<FieldMatrix<double,1,1>>;
  using CoefMatrices = CoefficientMatrices<MatrixType>;
public:
  GeometryBuffer() : bufferMap() {}

  std::pair<CoefMatrices&, bool>
  operator()(const GeometryComparable<Geometry>& geometry)
  {
    auto [pos, inserted] = bufferMap.emplace(geometry, CoefMatrices{});
    return {pos->second, !inserted};
  }

  unsigned int size() const
  {
    return bufferMap.size();
  }
private:
  std::map<GeometryComparable<Geometry>,
           CoefMatrices>                    bufferMap;
};



template<typename BilinForm, typename InnerProd>
class BufferedTestspaceCoefficientMatrix
{
  public:
  using BilinearForm = BilinForm;
  using InnerProduct = InnerProd;
  using SolutionSpaces = typename BilinearForm::SolutionSpaces;
  using TestSpaces     = typename BilinearForm::TestSpaces;
  using GridView = typename std::tuple_element_t<0, SolutionSpaces>::GridView;
  using Entity = typename GridView::template Codim<0>::Entity;
  using Geometry = decltype(std::declval<Entity>().geometry());
  using GeoBuffer = GeometryBuffer<Geometry>;

  private:
  using SolutionLocalViews = detail::getLocalViews_t<SolutionSpaces>;
  using TestLocalViews     = detail::getLocalViews_t<TestSpaces>;

  using MatrixType = Matrix<FieldMatrix<double,1,1>>;

  public:
  BufferedTestspaceCoefficientMatrix(BilinearForm& bilinForm,
                                     InnerProduct& innerProd,
                                     GeoBuffer& buffer) :
    bilinearForm_(bilinForm),
    innerProduct_(innerProd),
    localViewsSolution_(detail::getLocalViews(
                          *bilinearForm_.getSolutionSpaces())),
    localViewsTest_(detail::getLocalViews(*bilinearForm_.getTestSpaces())),
    geometryBuffer_(buffer)
  {}

  BufferedTestspaceCoefficientMatrix(BufferedTestspaceCoefficientMatrix&&)
    = default;

  BufferedTestspaceCoefficientMatrix(BufferedTestspaceCoefficientMatrix&)
    = delete;

  void bind(const Entity& e)
  {
    using namespace Dune::detail;
    GeometryComparable<Geometry> newGeometry;
    newGeometry.set(e.geometry());
    auto [coefMatrices, initialized] = geometryBuffer_(newGeometry);
    coefMatrices_ = &coefMatrices;
    if (!initialized)
    {
      coefMatrices_->set(e,
                         bilinearForm_,
                         innerProduct_,
                         localViewsSolution_,
                         localViewsTest_);
    }
  }

  const MatrixType& coefficientMatrix() const
  {
    return coefMatrices_->coefficientMatrix();
  }

  const MatrixType& systemMatrix() const
  {
    return coefMatrices_->systemMatrix();
  }

  const BilinearForm& bilinearForm() const
  {
    return bilinearForm_;
  }

  const InnerProduct& innerProduct() const
  {
    return innerProduct_;
  }

  const GridView& gridView() const
  {
    return std::get<0>(*bilinearForm_.getSolutionSpaces()).gridView();
  }

  private:
  BilinearForm&         bilinearForm_;
  InnerProduct&         innerProduct_;
  SolutionLocalViews    localViewsSolution_;
  TestLocalViews        localViewsTest_;
  CoefficientMatrices<MatrixType>* coefMatrices_;
  GeoBuffer&            geometryBuffer_;
};


} // end namespace Dune


#endif // DUNE_TESTSPACECOEFFICIENTMATRIX_HH
