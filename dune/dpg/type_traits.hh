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
  template<class size_type>
  class FlatMultiIndex;

  template<class NF>
  class DefaultGlobalBasis;

  template<typename GV, int s, int k, class MI, class ST>
  class PQkDGSubsampledDGNodeFactory;

  template<typename GV, int s, int k, class MI, class ST>
  class PQkSubsampledDGNodeFactory;

  template<typename GV, int k, class MI, class ST>
  class PQkTransportFactory;

  template<typename TestspaceCoefficientMatrix, std::size_t testIndex,
           class MI, class ST>
  class OptimalTestBasisNodeFactory;
}

template<class T, int size> class FieldVector;


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

template<class GV, int s, int k, class ST>
struct is_SubsampledFiniteElement<Functions::DefaultGlobalBasis<
             Functions::PQkDGSubsampledDGNodeFactory
                 <GV, s, k, Functions::FlatMultiIndex<ST>, ST> > >
       : std::true_type {};

template<class GV, int s, int k, class ST>
struct is_SubsampledFiniteElement<Functions::DefaultGlobalBasis<
             Functions::PQkSubsampledDGNodeFactory
                 <GV, s, k, Functions::FlatMultiIndex<ST>, ST> > >
       : std::true_type {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class ST>
struct is_SubsampledFiniteElement<Functions::DefaultGlobalBasis<
            Functions::OptimalTestBasisNodeFactory<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<ST>, ST> > >
  : is_SubsampledFiniteElement<typename std::tuple_element<testIndex,
                                 typename TestspaceCoefficientMatrix::TestSpaces
                               >::type> {};

template <typename FiniteElement>
struct numberOfSamples : std::integral_constant<int, 1> {};

template<class GV, int s, int k, class ST>
struct numberOfSamples<Functions::DefaultGlobalBasis<
             Functions::PQkDGSubsampledDGNodeFactory
                 <GV, s, k, Functions::FlatMultiIndex<ST>, ST> > >
       : std::integral_constant<int, s> {};

template<class GV, int s, int k, class ST>
struct numberOfSamples<Functions::DefaultGlobalBasis<
             Functions::PQkSubsampledDGNodeFactory
                 <GV, s, k, Functions::FlatMultiIndex<ST>, ST> > >
       : std::integral_constant<int, s> {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class ST>
struct numberOfSamples<Functions::DefaultGlobalBasis<
            Functions::OptimalTestBasisNodeFactory<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<ST>, ST> > >
  : numberOfSamples<typename std::tuple_element<testIndex,
                      typename TestspaceCoefficientMatrix::TestSpaces
                    >::type> {};

template<typename FiniteElement>
struct is_TransportFiniteElement : std::false_type {};

template<typename GV, int k, class ST>
struct is_TransportFiniteElement<Functions::DefaultGlobalBasis<
             Functions::PQkTransportFactory
                 <GV, k, Functions::FlatMultiIndex<ST>, ST> > >
           : std::true_type {};

template<typename TestspaceCoefficientMatrix, std::size_t testIndex, class ST>
struct is_TransportFiniteElement<Functions::DefaultGlobalBasis<
            Functions::OptimalTestBasisNodeFactory<
                TestspaceCoefficientMatrix, testIndex,
                Functions::FlatMultiIndex<ST>, ST> > >
  : is_TransportFiniteElement<typename std::tuple_element<testIndex,
                                typename TestspaceCoefficientMatrix::TestSpaces
                              >::type> {};

} // end namespace Dune

#endif // DUNE_DPG_TYPE_TRAITS_HH
