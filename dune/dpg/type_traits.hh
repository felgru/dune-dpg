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

  template<class NF>
  class DefaultGlobalBasis;

  template<typename GV, int k, class MI>
  class PQkNodeFactory;

  template<typename GV, int k, class MI>
  class LagrangeDGNodeFactory;

  template<typename GV, int level, int k, class MI>
  class PQkDGRefinedDGNodeFactory;

  template<typename GV, int s, int k, class MI>
  class PQkDGSubsampledDGNodeFactory;

  template<typename GV, int s, int k, class MI>
  class PQkSubsampledDGNodeFactory;

  template<typename GV, int k, class MI>
  class PQkTransportFactory;

  template<typename TestspaceCoefficientMatrix, std::size_t testIndex,
           class MI>
  class OptimalTestBasisNodeFactory;

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

/****************************
 * Traits for finite elements
 ****************************/

template <typename FiniteElement>
struct is_OptimalTestSpace : std::false_type {};

#ifndef DOXYGEN
template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct is_OptimalTestSpace<Functions::DefaultGlobalBasis<
            Functions::OptimalTestBasisNodeFactory<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> > > >
  : std::true_type {};
#endif

template <typename FiniteElement>
struct is_RefinedFiniteElement;

template <typename FiniteElement>
struct is_DGRefinedFiniteElement : std::false_type {};

#ifndef DOXYGEN
template<typename GV, int level, int k>
struct is_DGRefinedFiniteElement<Functions::DefaultGlobalBasis<
               Functions::PQkDGRefinedDGNodeFactory
                   <GV, level, k, Functions::FlatMultiIndex<std::size_t> > > >
       : std::true_type {};
#endif

template <typename FiniteElement>
struct is_ContinuouslyRefinedFiniteElement : std::false_type {};

#ifndef DOXYGEN
template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct is_ContinuouslyRefinedFiniteElement<Functions::DefaultGlobalBasis<
            Functions::OptimalTestBasisNodeFactory<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> > > >
  : is_RefinedFiniteElement<typename std::tuple_element<testIndex,
                              typename TestspaceCoefficientMatrix::TestSpaces
                            >::type> {};
#endif

template <typename FiniteElement>
struct is_RefinedFiniteElement
  : std::integral_constant<bool,
         is_DGRefinedFiniteElement<FiniteElement>::value
      || is_ContinuouslyRefinedFiniteElement<FiniteElement>::value> {};

template <typename FiniteElement>
struct levelOfFE : std::integral_constant<int, 0> {};

#ifndef DOXYGEN
template<class GV, int level, int k>
struct levelOfFE<Functions::DefaultGlobalBasis<
             Functions::PQkDGRefinedDGNodeFactory
                 <GV, level, k, Functions::FlatMultiIndex<std::size_t> > > >
       : std::integral_constant<int, level> {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct levelOfFE<Functions::DefaultGlobalBasis<
            Functions::OptimalTestBasisNodeFactory<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> > > >
  : levelOfFE<typename std::tuple_element<testIndex,
                typename TestspaceCoefficientMatrix::TestSpaces
              >::type> {};
#endif

template <typename FiniteElement>
struct is_SubsampledFiniteElement : std::false_type {};

#ifndef DOXYGEN
template<class GV, int s, int k>
struct is_SubsampledFiniteElement<Functions::DefaultGlobalBasis<
             Functions::PQkDGSubsampledDGNodeFactory
                 <GV, s, k, Functions::FlatMultiIndex<std::size_t> > > >
       : std::true_type {};

template<class GV, int s, int k>
struct is_SubsampledFiniteElement<Functions::DefaultGlobalBasis<
             Functions::PQkSubsampledDGNodeFactory
                 <GV, s, k, Functions::FlatMultiIndex<std::size_t> > > >
       : std::true_type {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct is_SubsampledFiniteElement<Functions::DefaultGlobalBasis<
            Functions::OptimalTestBasisNodeFactory<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> > > >
  : is_SubsampledFiniteElement<typename std::tuple_element<testIndex,
                                 typename TestspaceCoefficientMatrix::TestSpaces
                               >::type> {};
#endif

template <typename FiniteElement>
struct numberOfSamples : std::integral_constant<int, 1> {};

#ifndef DOXYGEN
template<class GV, int s, int k>
struct numberOfSamples<Functions::DefaultGlobalBasis<
             Functions::PQkDGSubsampledDGNodeFactory
                 <GV, s, k, Functions::FlatMultiIndex<std::size_t> > > >
       : std::integral_constant<int, s> {};

template<class GV, int s, int k>
struct numberOfSamples<Functions::DefaultGlobalBasis<
             Functions::PQkSubsampledDGNodeFactory
                 <GV, s, k, Functions::FlatMultiIndex<std::size_t> > > >
       : std::integral_constant<int, s> {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct numberOfSamples<Functions::DefaultGlobalBasis<
            Functions::OptimalTestBasisNodeFactory<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> > > >
  : numberOfSamples<typename std::tuple_element<testIndex,
                      typename TestspaceCoefficientMatrix::TestSpaces
                    >::type> {};
#endif

template<typename FiniteElement>
struct is_TransportFiniteElement : std::false_type {};

#ifndef DOXYGEN
template<typename GV, int k>
struct is_TransportFiniteElement<Functions::DefaultGlobalBasis<
             Functions::PQkTransportFactory
                 <GV, k, Functions::FlatMultiIndex<std::size_t> > > >
           : std::true_type {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct is_TransportFiniteElement<Functions::DefaultGlobalBasis<
            Functions::OptimalTestBasisNodeFactory<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> > > >
  : is_TransportFiniteElement<typename std::tuple_element<testIndex,
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

#ifndef DOXYGEN
template<class T, class GridView>
using changeGridView_t = typename changeGridView<T, GridView>::type;

template<class NF, class GridView>
struct changeGridView<Functions::DefaultGlobalBasis<NF>, GridView>
{
  typedef Functions::DefaultGlobalBasis<changeGridView_t<NF, GridView>> type;
};

template<typename GV, int k, class MI, class GridView>
struct changeGridView<Functions::PQkNodeFactory<GV, k, MI>, GridView>
{
  typedef Functions::PQkNodeFactory<GridView, k, MI> type;
};

template<typename GV, int k, class MI, class GridView>
struct changeGridView<Functions::LagrangeDGNodeFactory<GV, k, MI>, GridView>
{
  typedef Functions::LagrangeDGNodeFactory<GridView, k, MI> type;
};

template<typename GV, int level, int k, class MI, class GridView>
struct changeGridView<Functions::PQkDGRefinedDGNodeFactory<GV, level, k, MI>,
                      GridView>
{
  typedef Functions::PQkDGRefinedDGNodeFactory<GridView, level, k, MI> type;
};

template<typename GV, int s, int k, class MI, class GridView>
struct changeGridView<Functions::PQkDGSubsampledDGNodeFactory<GV, s, k, MI>,
                      GridView>
{
  typedef Functions::PQkDGSubsampledDGNodeFactory<GridView, s, k, MI> type;
};

template<typename GV, int s, int k, class MI, class GridView>
struct changeGridView<Functions::PQkSubsampledDGNodeFactory<GV, s, k, MI>,
                      GridView>
{
  typedef Functions::PQkSubsampledDGNodeFactory<GridView, s, k, MI> type;
};

template<typename GV, int k, class MI, class GridView>
struct changeGridView<Functions::PQkTransportFactory<GV, k, MI>, GridView>
{
  typedef Functions::PQkTransportFactory<GridView, k, MI> type;
};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex,
         class MI, class GridView>
struct changeGridView<Functions::OptimalTestBasisNodeFactory
                       <TestspaceCoefficientMatrix,
                        testIndex, MI>, GridView>
{
  typedef Functions::OptimalTestBasisNodeFactory<
    changeGridView_t<TestspaceCoefficientMatrix, GridView>,
    testIndex, MI>   type;
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



template<class T, class NewTestSpaces>
using replaceTestSpaces_t =
    std::decay_t<decltype(replaceTestSpaces(std::declval<T>(),
                          std::declval<NewTestSpaces>()))>;

} // end namespace Dune

#endif // DUNE_DPG_TYPE_TRAITS_HH
