// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_TESTSPACECOEFFICIENTMATRIX_HH
#define DUNE_TESTSPACECOEFFICIENTMATRIX_HH

#include <tuple>
#include <functional>
#include <memory>
#include <type_traits>

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

//#include <dune/localfunctions/optimaltestfunctions/optimaltest.hh>
//#include <dune/localfunctions/optimaltestfunctions/refinedoptimaltest.hh>

//#include <dune/typetree/leafnode.hh>

//#include <dune/functions/functionspacebases/nodes.hh>
//#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
//#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/cholesky.hh>
//#include <dune/dpg/type_traits.hh>


namespace Dune {

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
  }

  bool isSame(Geometry geometry)
  {
    // TODO: The return type of Geometry::corners() should be changed to unsigned.
    if (geometry.corners()==(int)corners)
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
  }

  private:
  unsigned int corners;
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
             >::type SolutionLocalViews;

  typedef typename boost::fusion::result_of::as_vector<
             typename boost::fusion::result_of::transform<
                         TestSpaces,
                         detail::getLocalView
                      >::type
             >::type TestLocalViews;

  typedef Matrix<FieldMatrix<double,1,1> > MatrixType;

  public:
  TestspaceCoefficientMatrix(BilinearForm& bilinForm, InnerProduct& innerProd) :
    bilinearForm_(bilinForm),
    innerProduct_(innerProd),
    gridView_(std::get<0>(bilinForm.getSolutionSpaces()).gridView()),
    localViewsSolution_(boost::fusion::as_vector(
                boost::fusion::transform(bilinearForm_.getSolutionSpaces(),
                                         detail::getLocalView()))),
    localViewsTest_(boost::fusion::as_vector(
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
      for_each(localViewsTest_, applyBind<decltype(e)>(e));
      for_each(localViewsSolution_, applyBind<decltype(e)>(e));

      MatrixType stiffnessMatrix;

      innerProduct_.bind(localViewsTest_);
      innerProduct_.getLocalMatrix(stiffnessMatrix);

      MatrixType bilinearMatrix;

      bilinearForm_.bind(localViewsTest_, localViewsSolution_);
      bilinearForm_.getLocalMatrix(bilinearMatrix);

      Cholesky<MatrixType> cholesky(stiffnessMatrix);
      cholesky.apply(bilinearMatrix, coefficientMatrix_);

      unsigned int m = coefficientMatrix_.M();
      localMatrix_.setSize(m, m);

      for (unsigned int i=0; i<m; i++)
      {
        for (unsigned int j=0; j<m; j++)
        {
          localMatrix_[i][j]=0;
          for (unsigned int k=0; k<coefficientMatrix_.N(); k++)
          {
            localMatrix_[i][j]+=(bilinearMatrix[k][i]*coefficientMatrix_[k][j]);
          }
        }
      }
      geometryBufferIsSet_ = true;
    }
  }

  MatrixType& coefficientMatrix()
  {
    return coefficientMatrix_;
  }

  MatrixType& localMatrix()
  {
    return localMatrix_;
  }

  BilinearForm& bilinearForm() const
  {
    return bilinearForm_;
  }

  InnerProduct& innerProduct() const
  {
    return innerProduct_;
  }

  const GridView& gridView() const
  {
    return gridView_;
  }

  private:
  BilinearForm&      bilinearForm_;
  InnerProduct&      innerProduct_;
  GridView           gridView_;
  SolutionLocalViews localViewsSolution_;
  TestLocalViews     localViewsTest_;
  GeometryBuffer<Geometry> geometryBuffer_;
  bool               geometryBufferIsSet_;
  MatrixType         coefficientMatrix_;        // G^{-1}B
  MatrixType         localMatrix_;              // B^TG^{-1}B
};



} // end namespace Dune


#endif // DUNE_TESTSPACECOEFFICIENTMATRIX_HH
