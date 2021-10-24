// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_TYPE_TRAITS_HH
#define DUNE_DPG_TYPE_TRAITS_HH

#include <type_traits>
#include <dune/common/keywords.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/version.hh>

#ifndef DOXYGEN
namespace std {
  template<class T, class Alloc> class vector;
  template< std::size_t I, class T > class tuple_element;
}
#endif

namespace Dune {

#ifndef DOXYGEN
namespace Functions {
#if DUNE_VERSION_LT(DUNE_FUNCTIONS,2,9)
  template<class size_type>
  class FlatMultiIndex;
#endif

  template<class PB>
  class DefaultGlobalBasis;

  template<class PB>
  class ConstrainedGlobalBasis;

  template<class PB>
  class RefinedGlobalBasis;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  template<typename GV, int k, typename R>
#else
  template<typename GV, int k, class MI, typename R>
#endif
  class LagrangePreBasis;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  template<typename GV, int k>
#else
  template<typename GV, int k, class MI>
#endif
  class LagrangeDGPreBasis;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  template<typename GV, int level, int k, typename R>
#else
  template<typename GV, int level, int k, class MI, typename R>
#endif
  class LagrangeDGRefinedDGPreBasis;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  template<typename GV, int k, typename R>
#else
  template<typename GV, int k, class MI, typename R>
#endif
  class BernsteinPreBasis;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  template<typename GV, int k, typename R>
#else
  template<typename GV, int k, class MI, typename R>
#endif
  class BernsteinDGPreBasis;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  template<typename GV, int level, int k, typename R>
#else
  template<typename GV, int level, int k, class MI, typename R>
#endif
  class BernsteinDGRefinedDGPreBasis;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  template<typename GV, int s, int k, typename R>
#else
  template<typename GV, int s, int k, class MI, typename R>
#endif
  class LagrangeDGSubsampledDGPreBasis;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  template<typename GV, int s, int k, typename R>
#else
  template<typename GV, int s, int k, class MI, typename R>
#endif
  class LagrangeSubsampledDGPreBasis;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  template<typename GV, int k, typename R>
#else
  template<typename GV, int k, class MI, typename R>
#endif
  class LagrangeTracePreBasis;

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
  template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
#else
  template<typename TestspaceCoefficientMatrix, std::size_t testIndex,
           class MI>
#endif
  class OptimalTestBasisPreBasis;

  template<typename GV, class MI, typename R>
  class HangingNodeBernsteinP2PreBasis;

  template<typename GV, class MI, typename R>
  class HangingNodeLagrangeP2PreBasis;

  template<typename InnerProduct>
  class NormalizedPreBasis;

  template<typename InnerProduct>
  class NormalizedRefinedPreBasis;

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
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex>
#else
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> >
#endif
            > >
  : std::true_type {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct is_OptimalTestSpace<Functions::RefinedGlobalBasis<
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex>
#else
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> >
#endif
            > >
  : std::true_type {};
#endif

template <typename GlobalBasis>
struct is_RefinedFiniteElement;

template <typename GlobalBasis>
struct is_DGRefinedFiniteElement : std::false_type {};

#ifndef DOXYGEN
template<typename GV, int level, int k, typename R>
struct is_DGRefinedFiniteElement<Functions::RefinedGlobalBasis<
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
               Functions::LagrangeDGRefinedDGPreBasis<GV, level, k, R>
#else
               Functions::LagrangeDGRefinedDGPreBasis
                   <GV, level, k, Functions::FlatMultiIndex<std::size_t>, R>
#endif
               >>
       : std::true_type {};

template<typename GV, int level, int k, typename R>
struct is_DGRefinedFiniteElement<Functions::RefinedGlobalBasis<
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
               Functions::BernsteinDGRefinedDGPreBasis<GV, level, k, R>
#else
               Functions::BernsteinDGRefinedDGPreBasis
                   <GV, level, k, Functions::FlatMultiIndex<std::size_t>, R>
#endif
               >>
       : std::true_type {};

template<typename InnerProduct>
struct is_DGRefinedFiniteElement<Functions::RefinedGlobalBasis<
                    Functions::NormalizedRefinedPreBasis<InnerProduct>>>
       : std::true_type {};
#endif

template <typename GlobalBasis>
struct is_ContinuouslyRefinedFiniteElement : std::false_type {};

#ifndef DOXYGEN
template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct is_ContinuouslyRefinedFiniteElement<Functions::RefinedGlobalBasis<
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex>
#else
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> >
#endif
            > >
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
template<class GV, int level, int k, typename R>
struct levelOfFE<Functions::RefinedGlobalBasis<
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
             Functions::LagrangeDGRefinedDGPreBasis<GV, level, k, R>
#else
             Functions::LagrangeDGRefinedDGPreBasis
                 <GV, level, k, Functions::FlatMultiIndex<std::size_t>, R>
#endif
             >>
       : std::integral_constant<int, level> {};

template<class GV, int level, int k, typename R>
struct levelOfFE<Functions::RefinedGlobalBasis<
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
             Functions::BernsteinDGRefinedDGPreBasis<GV, level, k, R>
#else
             Functions::BernsteinDGRefinedDGPreBasis
                 <GV, level, k, Functions::FlatMultiIndex<std::size_t>, R>
#endif
             >>
       : std::integral_constant<int, level> {};

template<class InnerProduct>
struct levelOfFE<Functions::RefinedGlobalBasis<
             Functions::NormalizedRefinedPreBasis<InnerProduct>>>
       : levelOfFE<typename std::tuple_element<static_cast<size_t>(0),
                              typename InnerProduct::TestSpaces>::type> {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct levelOfFE<Functions::RefinedGlobalBasis<
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex>
#else
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> >
#endif
            > >
  : levelOfFE<typename std::tuple_element<testIndex,
                typename TestspaceCoefficientMatrix::TestSpaces
              >::type> {};
#endif

template <typename GlobalBasis>
struct is_SubsampledFiniteElement : std::false_type {};

#ifndef DOXYGEN
template<class GV, int s, int k, typename R>
struct is_SubsampledFiniteElement<Functions::DefaultGlobalBasis<
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
             Functions::LagrangeDGSubsampledDGPreBasis<GV, s, k, R>
#else
             Functions::LagrangeDGSubsampledDGPreBasis
                 <GV, s, k, Functions::FlatMultiIndex<std::size_t>, R>
#endif
             >>
       : std::true_type {};

template<class GV, int s, int k, typename R>
struct is_SubsampledFiniteElement<Functions::DefaultGlobalBasis<
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
             Functions::LagrangeSubsampledDGPreBasis<GV, s, k, R>
#else
             Functions::LagrangeSubsampledDGPreBasis
                 <GV, s, k, Functions::FlatMultiIndex<std::size_t>, R>
#endif
             >>
       : std::true_type {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct is_SubsampledFiniteElement<Functions::DefaultGlobalBasis<
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex>
#else
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> >
#endif
            > >
  : is_SubsampledFiniteElement<typename std::tuple_element<testIndex,
                                 typename TestspaceCoefficientMatrix::TestSpaces
                               >::type> {};
#endif

template <typename GlobalBasis>
struct numberOfSamples : std::integral_constant<int, 1> {};

#ifndef DOXYGEN
template<class GV, int s, int k, typename R>
struct numberOfSamples<Functions::DefaultGlobalBasis<
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
             Functions::LagrangeDGSubsampledDGPreBasis<GV, s, k, R>
#else
             Functions::LagrangeDGSubsampledDGPreBasis
                 <GV, s, k, Functions::FlatMultiIndex<std::size_t>, R>
#endif
             >>
       : std::integral_constant<int, s> {};

template<class GV, int s, int k, typename R>
struct numberOfSamples<Functions::DefaultGlobalBasis<
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
             Functions::LagrangeSubsampledDGPreBasis<GV, s, k, R>
#else
             Functions::LagrangeSubsampledDGPreBasis
                 <GV, s, k, Functions::FlatMultiIndex<std::size_t>, R>
#endif
             >>
       : std::integral_constant<int, s> {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct numberOfSamples<Functions::DefaultGlobalBasis<
#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex>
#else
            Functions::OptimalTestBasisPreBasis<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<std::size_t> >
#endif
            > >
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

template<class PB, class GridView>
struct changeGridView<Functions::RefinedGlobalBasis<PB>, GridView>
{
  typedef Functions::RefinedGlobalBasis<changeGridView_t<PB, GridView>> type;
};

#if DUNE_VERSION_GTE(DUNE_FUNCTIONS,2,9)
template<typename GV, int k, typename R, class GridView>
struct changeGridView<Functions::LagrangePreBasis<GV, k, R>, GridView>
{
  typedef Functions::LagrangePreBasis<GridView, k, R> type;
};

template<typename GV, int k, class GridView>
struct changeGridView<Functions::LagrangeDGPreBasis<GV, k>, GridView>
{
  typedef Functions::LagrangeDGPreBasis<GridView, k> type;
};

template<typename GV, int level, int k, typename R, class GridView>
struct changeGridView<Functions::LagrangeDGRefinedDGPreBasis<GV, level, k, R>,
                      GridView>
{
  typedef Functions::LagrangeDGRefinedDGPreBasis<GridView, level, k, R> type;
};

template<typename GV, int k, typename R, class GridView>
struct changeGridView<Functions::BernsteinPreBasis<GV, k, R>, GridView>
{
  typedef Functions::BernsteinPreBasis<GridView, k, R> type;
};

template<typename GV, int k, typename R, class GridView>
struct changeGridView<Functions::BernsteinDGPreBasis<GV, k, R>, GridView>
{
  typedef Functions::BernsteinDGPreBasis<GridView, k, R> type;
};

template<typename GV, int level, int k, typename R, class GridView>
struct changeGridView<Functions::BernsteinDGRefinedDGPreBasis
                                                <GV, level, k, R>,
                      GridView>
{
  typedef Functions::BernsteinDGRefinedDGPreBasis<GridView, level, k, R>
      type;
};

template<typename GV, int s, int k, typename R, class GridView>
struct changeGridView<Functions::LagrangeDGSubsampledDGPreBasis<GV, s, k, R>,
                      GridView>
{
  typedef Functions::LagrangeDGSubsampledDGPreBasis<GridView, s, k, R> type;
};

template<typename GV, int s, int k, typename R, class GridView>
struct changeGridView<Functions::LagrangeSubsampledDGPreBasis<GV, s, k, R>,
                      GridView>
{
  typedef Functions::LagrangeSubsampledDGPreBasis<GridView, s, k, R> type;
};

template<typename GV, int k, typename R, class GridView>
struct changeGridView<Functions::LagrangeTracePreBasis<GV, k, R>,
                      GridView>
{
  typedef Functions::LagrangeTracePreBasis<GridView, k, R> type;
};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex,
         class GridView>
struct changeGridView<Functions::OptimalTestBasisPreBasis
                       <TestspaceCoefficientMatrix,
                        testIndex>, GridView>
{
  typedef Functions::OptimalTestBasisPreBasis<
    changeGridView_t<TestspaceCoefficientMatrix, GridView>,
    testIndex>   type;
};
#else
template<typename GV, int k, class MI, typename R, class GridView>
struct changeGridView<Functions::LagrangePreBasis<GV, k, MI, R>, GridView>
{
  typedef Functions::LagrangePreBasis<GridView, k, MI, R> type;
};

template<typename GV, int k, class MI, class GridView>
struct changeGridView<Functions::LagrangeDGPreBasis<GV, k, MI>, GridView>
{
  typedef Functions::LagrangeDGPreBasis<GridView, k, MI> type;
};

template<typename GV, int level, int k, class MI, typename R, class GridView>
struct changeGridView<Functions::LagrangeDGRefinedDGPreBasis<GV, level, k, MI, R>,
                      GridView>
{
  typedef Functions::LagrangeDGRefinedDGPreBasis<GridView, level, k, MI, R> type;
};

template<typename GV, int k, class MI, typename R, class GridView>
struct changeGridView<Functions::BernsteinPreBasis<GV, k, MI, R>, GridView>
{
  typedef Functions::BernsteinPreBasis<GridView, k, MI, R> type;
};

template<typename GV, int k, class MI, typename R, class GridView>
struct changeGridView<Functions::BernsteinDGPreBasis<GV, k, MI, R>, GridView>
{
  typedef Functions::BernsteinDGPreBasis<GridView, k, MI, R> type;
};

template<typename GV, int level, int k, class MI, typename R, class GridView>
struct changeGridView<Functions::BernsteinDGRefinedDGPreBasis
                                                <GV, level, k, MI, R>,
                      GridView>
{
  typedef Functions::BernsteinDGRefinedDGPreBasis<GridView, level, k, MI, R>
      type;
};

template<typename GV, int s, int k, class MI, typename R, class GridView>
struct changeGridView<Functions::LagrangeDGSubsampledDGPreBasis<GV, s, k, MI, R>,
                      GridView>
{
  typedef Functions::LagrangeDGSubsampledDGPreBasis<GridView, s, k, MI, R> type;
};

template<typename GV, int s, int k, class MI, typename R, class GridView>
struct changeGridView<Functions::LagrangeSubsampledDGPreBasis<GV, s, k, MI, R>,
                      GridView>
{
  typedef Functions::LagrangeSubsampledDGPreBasis<GridView, s, k, MI, R> type;
};

template<typename GV, int k, class MI, typename R, class GridView>
struct changeGridView<Functions::LagrangeTracePreBasis<GV, k, MI, R>,
                      GridView>
{
  typedef Functions::LagrangeTracePreBasis<GridView, k, MI, R> type;
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
#endif

template<typename GV, class MI, typename R, class GridView>
struct changeGridView<Functions::HangingNodeBernsteinP2PreBasis<GV, MI, R>,
                      GridView>
{
  typedef Functions::HangingNodeBernsteinP2PreBasis<GridView, MI, R> type;
};

template<typename GV, class MI, typename R, class GridView>
struct changeGridView<Functions::HangingNodeLagrangeP2PreBasis<GV, MI, R>,
                      GridView>
{
  typedef Functions::HangingNodeLagrangeP2PreBasis<GridView, MI, R> type;
};

template<typename InnerProduct, class GridView>
struct changeGridView<Functions::NormalizedPreBasis<InnerProduct>,
                      GridView>
{
  typedef Functions::NormalizedPreBasis<
            changeGridView_t<InnerProduct, GridView>> type;
};

template<typename InnerProduct, class GridView>
struct changeGridView<Functions::NormalizedRefinedPreBasis<InnerProduct>,
                      GridView>
{
  typedef Functions::NormalizedRefinedPreBasis<
            changeGridView_t<InnerProduct, GridView>> type;
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


template<class Term>
struct uses_only_constant_coefficients : std::false_type {};

template<class... Terms>
struct uses_only_constant_coefficients<std::tuple<Terms...>>
  : std::conjunction<
      uses_only_constant_coefficients<typename Terms::Term>...> {};

template<class Term>
DUNE_INLINE_VARIABLE constexpr bool uses_only_constant_coefficients_v
    = uses_only_constant_coefficients<Term>::value;

} // end namespace Dune

#endif // DUNE_DPG_TYPE_TRAITS_HH
