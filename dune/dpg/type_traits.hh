// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_TYPE_TRAITS_HH
#define DUNE_DPG_TYPE_TRAITS_HH

#include<type_traits>

namespace std {
  template<class T, class Alloc> class vector;
}

namespace Dune {

namespace Functions {
  template<class GV, int s, int k>
  class PQkSubsampledDGBasis;

  template<typename GV, int k>
  class PQkTransportBasis;

  template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
  class OptimalTestBasis;
}

template<class T, int size> class FieldVector;

/* Convenience interface for enable_if, taken from
 * http://flamingdangerzone.com/cxx11/2012/06/01/almost-static-if.html */

namespace detail {
  enum class enabler {};
}

template <typename Condition>
using EnableIf =
  typename std::enable_if<Condition::value, detail::enabler>::type;


template<class T>
struct is_vector : std::false_type {};

template<class T, class Alloc>
struct is_vector<std::vector<T, Alloc>> : std::true_type {};

template<class T, int size>
struct is_vector<Dune::FieldVector<T, size>> : std::true_type {};

/****************************
 * Traits for finite elements
 ****************************/

template <typename FiniteElement>
struct is_SubsampledFiniteElement : std::false_type {};

template<class GV, int s, int k>
struct is_SubsampledFiniteElement<Functions::PQkSubsampledDGBasis
                                  <GV, s, k> > : std::true_type {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct is_SubsampledFiniteElement<Functions::OptimalTestBasis
                                  <TestspaceCoefficientMatrix, testIndex> >
  : std::integral_constant<bool,
          is_SubsampledFiniteElement<typename std::tuple_element<testIndex,
                                     typename Functions::OptimalTestBasis
                                     <TestspaceCoefficientMatrix, testIndex>
                                     ::EnrichedTestspaces>::type>::value> {};

template <typename FiniteElement>
struct numberOfSamples : std::integral_constant<int, 1> {};

template<class GV, int s, int k>
struct numberOfSamples<Functions::PQkSubsampledDGBasis<GV, s, k> >
                : std::integral_constant<int, s> {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct numberOfSamples<Functions::OptimalTestBasis
                                  <TestspaceCoefficientMatrix, testIndex> >
                : std::integral_constant<int, numberOfSamples<
                                     typename std::tuple_element<testIndex,
                                     typename Functions::OptimalTestBasis
                                     <TestspaceCoefficientMatrix, testIndex>
                                     ::EnrichedTestspaces>::type>::value> {};

template<typename FiniteElement>
struct is_TransportFiniteElement : std::false_type {};

template<typename GV, int k>
struct is_TransportFiniteElement<Functions::PQkTransportBasis<GV, k> >
           : std::true_type {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex>
struct is_TransportFiniteElement<Functions::OptimalTestBasis
                                  <TestspaceCoefficientMatrix, testIndex> >
  : std::integral_constant<bool,
          is_TransportFiniteElement<typename std::tuple_element<testIndex,
                                    typename Functions::OptimalTestBasis
                                    <TestspaceCoefficientMatrix, testIndex>
                                    ::EnrichedTestspaces>::type>::value> {};

} // end namespace Dune

#endif // DUNE_DPG_TYPE_TRAITS_HH
