// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_TYPE_TRAITS_HH
#define DUNE_DPG_TYPE_TRAITS_HH

#include <type_traits>
#include <dune/common/tupleutility.hh>

#ifndef DOXYGEN
namespace std {
  template<class T, class Alloc> class vector;
  template< std::size_t I, class T > class tuple_element;
}
#endif

namespace Dune {

#ifndef DOXYGEN
namespace Functions {
  template<class size_type>
  class FlatMultiIndex;

  template<class PB>
  class DefaultGlobalBasis;

  template<class PB>
  class ConstrainedGlobalBasis;

  template<typename GV, int k, class MI>
  class PQkPreBasis;

  template<typename GV, int k, class MI>
  class LagrangeDGPreBasis;

  template<typename GV, int level, int k, class MI>
  class PQkDGRefinedDGPreBasis;

  template<typename GV, int k, class MI>
  class BernsteinPreBasis;

  template<typename GV, int k, class MI>
  class BernsteinDGPreBasis;

  template<typename GV, int level, int k, class MI>
  class BernsteinDGRefinedDGPreBasis;

  template<typename GV, int s, int k, class MI>
  class PQkDGSubsampledDGPreBasis;

  template<typename GV, int s, int k, class MI>
  class PQkSubsampledDGPreBasis;

  template<typename GV, int k, class MI>
  class PQkTracePreBasis;

  template<typename TestspaceCoefficientMatrix, std::size_t testIndex,
           class MI>
  class OptimalTestBasisPreBasis;

  template<typename GV, class MI>
  class HangingNodeP2PreBasis;

  template<typename BilinForm, typename InnerProd>
  class TestspaceCoefficientMatrix;
}

template<class TSpaces, class SolSpaces, class BilinearTerms>
class BilinearForm;

template<class TSpaces, class InnerProductTerms>
class InnerProduct;

template<class T, int size> class FieldVector;
#endif


#ifdef DOXYGEN
//! check if T is a std or Dune vector type
//!
//! is_vector<T>::value is true if T is a vector type and false otherwise.
template<class T>
struct is_vector {};
#else
template<class T>
struct is_vector : std::false_type {};

template<class T, class Alloc>
struct is_vector<std::vector<T, Alloc>> : std::true_type {};

template<class T, int size>
struct is_vector<Dune::FieldVector<T, size>> : std::true_type {};
#endif

template<class T>
struct tuple_size {};

template<class T, int size>
struct tuple_size<Dune::FieldVector<T, size>>
: std::integral_constant<size_t, size> {};

/****************************
 * Traits for finite elements
 ****************************/

template <typename GlobalBasis>
struct is_OptimalTestSpace : std::false_type {};

#ifndef DOXYGEN
template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct is_OptimalTestSpace<Functions::DefaultGlobalBasis<
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> > > >
  : std::true_type {};
#endif

template <typename GlobalBasis>
struct is_RefinedFiniteElement;

template <typename GlobalBasis>
struct is_DGRefinedFiniteElement : std::false_type {};

#ifndef DOXYGEN
template<typename GV, int level, int k>
struct is_DGRefinedFiniteElement<Functions::DefaultGlobalBasis<
               Functions::PQkDGRefinedDGPreBasis
                   <GV, level, k, Functions::FlatMultiIndex<std::size_t> > > >
       : std::true_type {};

template<typename GV, int level, int k>
struct is_DGRefinedFiniteElement<Functions::DefaultGlobalBasis<
               Functions::BernsteinDGRefinedDGPreBasis
                   <GV, level, k, Functions::FlatMultiIndex<std::size_t> > > >
       : std::true_type {};
#endif

template <typename GlobalBasis>
struct is_ContinuouslyRefinedFiniteElement : std::false_type {};

#ifndef DOXYGEN
template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct is_ContinuouslyRefinedFiniteElement<Functions::DefaultGlobalBasis<
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> > > >
  : is_RefinedFiniteElement<typename std::tuple_element<testIndex,
                              typename TestspaceCoefficientMatrix::TestSpaces
                            >::type> {};
#endif

template <typename GlobalBasis>
struct is_RefinedFiniteElement
  : std::integral_constant<bool,
         is_DGRefinedFiniteElement<GlobalBasis>::value
      || is_ContinuouslyRefinedFiniteElement<GlobalBasis>::value> {};

template <typename GlobalBasis>
struct levelOfFE : std::integral_constant<int, 0> {};

#ifndef DOXYGEN
template<class GV, int level, int k>
struct levelOfFE<Functions::DefaultGlobalBasis<
             Functions::PQkDGRefinedDGPreBasis
                 <GV, level, k, Functions::FlatMultiIndex<std::size_t> > > >
       : std::integral_constant<int, level> {};

template<class GV, int level, int k>
struct levelOfFE<Functions::DefaultGlobalBasis<
             Functions::BernsteinDGRefinedDGPreBasis
                 <GV, level, k, Functions::FlatMultiIndex<std::size_t> > > >
       : std::integral_constant<int, level> {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct levelOfFE<Functions::DefaultGlobalBasis<
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> > > >
  : levelOfFE<typename std::tuple_element<testIndex,
                typename TestspaceCoefficientMatrix::TestSpaces
              >::type> {};
#endif

template <typename GlobalBasis>
struct is_SubsampledFiniteElement : std::false_type {};

#ifndef DOXYGEN
template<class GV, int s, int k>
struct is_SubsampledFiniteElement<Functions::DefaultGlobalBasis<
             Functions::PQkDGSubsampledDGPreBasis
                 <GV, s, k, Functions::FlatMultiIndex<std::size_t> > > >
       : std::true_type {};

template<class GV, int s, int k>
struct is_SubsampledFiniteElement<Functions::DefaultGlobalBasis<
             Functions::PQkSubsampledDGPreBasis
                 <GV, s, k, Functions::FlatMultiIndex<std::size_t> > > >
       : std::true_type {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct is_SubsampledFiniteElement<Functions::DefaultGlobalBasis<
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> > > >
  : is_SubsampledFiniteElement<typename std::tuple_element<testIndex,
                                 typename TestspaceCoefficientMatrix::TestSpaces
                               >::type> {};
#endif

template <typename GlobalBasis>
struct numberOfSamples : std::integral_constant<int, 1> {};

#ifndef DOXYGEN
template<class GV, int s, int k>
struct numberOfSamples<Functions::DefaultGlobalBasis<
             Functions::PQkDGSubsampledDGPreBasis
                 <GV, s, k, Functions::FlatMultiIndex<std::size_t> > > >
       : std::integral_constant<int, s> {};

template<class GV, int s, int k>
struct numberOfSamples<Functions::DefaultGlobalBasis<
             Functions::PQkSubsampledDGPreBasis
                 <GV, s, k, Functions::FlatMultiIndex<std::size_t> > > >
       : std::integral_constant<int, s> {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct numberOfSamples<Functions::DefaultGlobalBasis<
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> > > >
  : numberOfSamples<typename std::tuple_element<testIndex,
                      typename TestspaceCoefficientMatrix::TestSpaces
                    >::type> {};
#endif

/*****************************************
 * Change the GridView of a nodal basis
 *****************************************/

//! change the GridView of a GlobalBasis, BilinearForm, etc.
//!
//! changeGridView<T, GV>::type changes the GridView template parameter
//! of T to GV.
template<class T, class GridView>
struct changeGridView {};

template<class T, class GridView>
using changeGridView_t = typename changeGridView<T, GridView>::type;

#ifndef DOXYGEN
template<class PB, class GridView>
struct changeGridView<Functions::DefaultGlobalBasis<PB>, GridView>
{
  typedef Functions::DefaultGlobalBasis<changeGridView_t<PB, GridView>> type;
};

template<class PB, class GridView>
struct changeGridView<Functions::ConstrainedGlobalBasis<PB>, GridView>
{
  typedef Functions::ConstrainedGlobalBasis<changeGridView_t<PB, GridView>>
        type;
};

template<typename GV, int k, class MI, class GridView>
struct changeGridView<Functions::PQkPreBasis<GV, k, MI>, GridView>
{
  typedef Functions::PQkPreBasis<GridView, k, MI> type;
};

template<typename GV, int k, class MI, class GridView>
struct changeGridView<Functions::LagrangeDGPreBasis<GV, k, MI>, GridView>
{
  typedef Functions::LagrangeDGPreBasis<GridView, k, MI> type;
};

template<typename GV, int level, int k, class MI, class GridView>
struct changeGridView<Functions::PQkDGRefinedDGPreBasis<GV, level, k, MI>,
                      GridView>
{
  typedef Functions::PQkDGRefinedDGPreBasis<GridView, level, k, MI> type;
};

template<typename GV, int k, class MI, class GridView>
struct changeGridView<Functions::BernsteinPreBasis<GV, k, MI>, GridView>
{
  typedef Functions::BernsteinPreBasis<GridView, k, MI> type;
};

template<typename GV, int k, class MI, class GridView>
struct changeGridView<Functions::BernsteinDGPreBasis<GV, k, MI>, GridView>
{
  typedef Functions::BernsteinDGPreBasis<GridView, k, MI> type;
};

template<typename GV, int level, int k, class MI, class GridView>
struct changeGridView<Functions::BernsteinDGRefinedDGPreBasis
                                                <GV, level, k, MI>,
                      GridView>
{
  typedef Functions::BernsteinDGRefinedDGPreBasis<GridView, level, k, MI>
      type;
};

template<typename GV, int s, int k, class MI, class GridView>
struct changeGridView<Functions::PQkDGSubsampledDGPreBasis<GV, s, k, MI>,
                      GridView>
{
  typedef Functions::PQkDGSubsampledDGPreBasis<GridView, s, k, MI> type;
};

template<typename GV, int s, int k, class MI, class GridView>
struct changeGridView<Functions::PQkSubsampledDGPreBasis<GV, s, k, MI>,
                      GridView>
{
  typedef Functions::PQkSubsampledDGPreBasis<GridView, s, k, MI> type;
};

template<typename GV, int k, class MI, class GridView>
struct changeGridView<Functions::PQkTracePreBasis<GV, k, MI>,
                      GridView>
{
  typedef Functions::PQkTracePreBasis<GridView, k, MI> type;
};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex,
         class MI, class GridView>
struct changeGridView<Functions::OptimalTestBasisPreBasis
                       <TestspaceCoefficientMatrix,
                        testIndex, MI>, GridView>
{
  typedef Functions::OptimalTestBasisPreBasis<
    changeGridView_t<TestspaceCoefficientMatrix, GridView>,
    testIndex, MI>   type;
};

template<typename GV, class MI, class GridView>
struct changeGridView<Functions::HangingNodeP2PreBasis<GV, MI>, GridView>
{
  typedef Functions::HangingNodeP2PreBasis<GridView, MI> type;
};

template<typename BilinForm, typename InnerProd, class GridView>
struct changeGridView<Functions::TestspaceCoefficientMatrix<
                          BilinForm, InnerProd>, GridView>
{
  typedef Functions::TestspaceCoefficientMatrix
      < changeGridView_t<BilinForm, GridView>
      , changeGridView_t<InnerProd, GridView>
      >  type;
};

template<class GridView, class... Spaces>
struct changeGridView<std::tuple<Spaces...>, GridView>
{
  template<class Space>
  struct changeGridView_
  { typedef changeGridView_t<Space, GridView> Type; };

  typedef typename ForEachType<changeGridView_, std::tuple<Spaces...>>::Type
      type;
};

template<class TSpaces, class SolSpaces, class BilinearTerms, class GridView>
struct changeGridView<
    BilinearForm<TSpaces, SolSpaces, BilinearTerms>,
    GridView>
{
  typedef BilinearForm
    < changeGridView_t<TSpaces, GridView>
    , changeGridView_t<SolSpaces, GridView>
    , BilinearTerms
    >  type;
};

template<class TSpaces, class InnerProductTerms, class GridView>
struct changeGridView<InnerProduct<TSpaces, InnerProductTerms>, GridView>
{
  typedef InnerProduct<changeGridView_t<TSpaces, GridView>, InnerProductTerms>
      type;
};
#endif

} // end namespace Dune

#endif // DUNE_DPG_TYPE_TRAITS_HH
