// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_BILINEARFORMFACTORY_HH
#define DUNE_DPG_BILINEARFORMFACTORY_HH

#include <tuple>
#include <type_traits>

#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/bilinearform.hh>
#include <dune/dpg/integralterm.hh>
#include <dune/dpg/spacetuple.hh>

namespace Dune {
  template<class TSpaces, class SolSpaces, class... BilinearTerms>
  class BilinearFormFactory;

  template<class TSpaces, class SolSpaces>
  BilinearFormFactory<TSpaces, SolSpaces>
  bilinearFormWithSpaces(const TSpaces& testSpaces,
                         const SolSpaces& solutionSpaces);

  template<class TSpaces, class SolSpaces, class... BilinearTerms>
  class BilinearFormFactory {
    static_assert(is_SpaceTuplePtr<TSpaces>::value,
        "TSpaces needs to be a SpaceTuplePtr!");
    static_assert(is_SpaceTuplePtr<SolSpaces>::value,
        "SolSpaces needs to be a SpaceTuplePtr!");
    using TestSpacesPtr = TSpaces;
    using SolutionSpacesPtr = SolSpaces;
    using TestSpaces = typename TestSpacesPtr::element_type;
    using SolutionSpaces = typename SolutionSpacesPtr::element_type;
    using BilinearTermsTuple = std::tuple<BilinearTerms...>;

    BilinearFormFactory(const TestSpacesPtr& testSpaces,
                        const SolutionSpacesPtr& solutionSpaces,
                        const BilinearTermsTuple& terms)
      : testSpaces(testSpaces)
      , solutionSpaces(solutionSpaces)
      , terms(terms)
      {}

    friend
    BilinearFormFactory<TSpaces, SolSpaces>
    bilinearFormWithSpaces<TSpaces, SolSpaces>(const TSpaces& testSpaces,
                           const SolSpaces& solutionSpaces);

    template<class TSpaces_, class SolSpaces_, class... BilinearTerms_>
    friend class BilinearFormFactory;

    TestSpacesPtr testSpaces;
    SolutionSpacesPtr solutionSpaces;
    BilinearTermsTuple terms;

  public:
    template<size_t lhsSpaceIndex,
             size_t rhsSpaceIndex,
             IntegrationType integrationType,
             DomainOfIntegration domainOfIntegration,
             class Factor,
             typename std::enable_if<
                         integrationType == IntegrationType::valueValue
                      || integrationType == IntegrationType::normalSign>::type*
                    = nullptr
            >
    auto addIntegralTerm(Factor c)
    {
      auto newTerm = make_IntegralTerm<lhsSpaceIndex,
                                       rhsSpaceIndex,
                                       integrationType,
                                       domainOfIntegration>(c);
      using NewTerm = decltype(newTerm);
      using NewBilinearFormFactory
        = BilinearFormFactory<TSpaces, SolSpaces, BilinearTerms..., NewTerm>;
      return NewBilinearFormFactory{
               testSpaces, solutionSpaces,
               std::tuple_cat(terms, std::make_tuple(std::move(newTerm)))
             };
    }

    template<size_t lhsSpaceIndex,
             size_t rhsSpaceIndex,
             IntegrationType integrationType,
             DomainOfIntegration domainOfIntegration,
             class Factor, class Direction,
             typename std::enable_if<
                         integrationType == IntegrationType::gradValue
                      || integrationType == IntegrationType::valueGrad
                      || integrationType == IntegrationType::gradGrad
                      || integrationType == IntegrationType::normalVector
                      || integrationType ==
                                IntegrationType::travelDistanceWeighted
                      >::type*
               = nullptr
            >
    auto addIntegralTerm(Factor c, Direction beta)
    {
      auto newTerm = make_IntegralTerm<lhsSpaceIndex,
                                       rhsSpaceIndex,
                                       integrationType,
                                       domainOfIntegration>(c, beta);
      using NewTerm = decltype(newTerm);
      using NewBilinearFormFactory
        = BilinearFormFactory<TSpaces, SolSpaces, BilinearTerms..., NewTerm>;
      return NewBilinearFormFactory{
               testSpaces, solutionSpaces,
               std::tuple_cat(terms, std::make_tuple(std::move(newTerm)))
             };
    }

    template<size_t lhsSpaceIndex,
             size_t rhsSpaceIndex,
             IntegrationType integrationType,
             DomainOfIntegration domainOfIntegration,
             class Factor, class Direction,
             typename std::enable_if<
                         integrationType == IntegrationType::gradGrad>::type*
                    = nullptr
            >
    auto addIntegralTerm(Factor c, Direction lhsBeta, Direction rhsBeta)
    {
      auto newTerm = make_IntegralTerm<lhsSpaceIndex,
                                       rhsSpaceIndex,
                                       integrationType,
                                       domainOfIntegration>
                                      (c, lhsBeta, rhsBeta);
      using NewTerm = decltype(newTerm);
      using NewBilinearFormFactory
        = BilinearFormFactory<TSpaces, SolSpaces, BilinearTerms..., NewTerm>;
      return NewBilinearFormFactory{
               testSpaces, solutionSpaces,
               std::tuple_cat(terms, std::make_tuple(std::move(newTerm)))
             };
    }

    BilinearForm<TSpaces, SolSpaces, std::tuple<BilinearTerms...>>
    create()
    {
      return {testSpaces, solutionSpaces, terms};
    }
  };

  template<class TSpaces, class SolSpaces>
  BilinearFormFactory<TSpaces, SolSpaces>
  bilinearFormWithSpaces(const TSpaces& testSpaces,
                         const SolSpaces& solutionSpaces)
  {
    return {testSpaces, solutionSpaces, std::tuple<>{}};
  }
}

#endif
