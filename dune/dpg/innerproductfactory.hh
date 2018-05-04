// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_INNERPRODUCTFACTORY_HH
#define DUNE_DPG_INNERPRODUCTFACTORY_HH

#include <tuple>
#include <type_traits>

#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/functions/gridviewfunctions.hh>
#include <dune/dpg/innerproduct.hh>
#include <dune/dpg/integralterm.hh>
#include <dune/dpg/spacetuple.hh>

namespace Dune {
  template<class TSpaces, class... InnerProductTerms>
  class InnerProductFactory;

  template<class TSpaces>
  InnerProductFactory<TSpaces>
  innerProductWithSpace(const TSpaces& testSpaces);

  template<class TSpaces, class... InnerProductTerms>
  class InnerProductFactory {
    static_assert(is_SpaceTuplePtr<TSpaces>::value,
        "TSpaces needs to be a SpaceTuplePtr!");
    using TestSpacesPtr = TSpaces;
    using TestSpaces = typename TestSpacesPtr::element_type;
    using InnerProductTermsTuple = std::tuple<InnerProductTerms...>;
    using GridView = typename std::tuple_element_t<0,TestSpaces>::GridView;

    InnerProductFactory(const TestSpacesPtr& testSpaces,
                        const InnerProductTermsTuple& terms)
      : testSpaces(testSpaces)
      , terms(terms)
      {}

    friend
    InnerProductFactory<TSpaces>
    innerProductWithSpace<TSpaces>(const TSpaces& testSpaces);

    template<class TSpaces_, class... InnerProductTerms_>
    friend class InnerProductFactory;

    TestSpacesPtr testSpaces;
    InnerProductTermsTuple terms;

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
      auto cFunc = Functions::detail::toGridViewFunction<GridView>(c);
      auto newTerm = make_IntegralTerm<lhsSpaceIndex,
                                       rhsSpaceIndex,
                                       integrationType,
                                       domainOfIntegration>(std::move(cFunc));
      using NewTerm = decltype(newTerm);
      using NewInnerProductFactory
        = InnerProductFactory<TSpaces, InnerProductTerms..., NewTerm>;
      return NewInnerProductFactory{
               testSpaces,
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
      auto cFunc = Functions::detail::toGridViewFunction<GridView>(c);
      auto betaFunc = Functions::detail::toGridViewFunction<GridView>(beta);
      auto newTerm = make_IntegralTerm<lhsSpaceIndex,
                                       rhsSpaceIndex,
                                       integrationType,
                                       domainOfIntegration>
                                      (std::move(cFunc), std::move(betaFunc));
      using NewTerm = decltype(newTerm);
      using NewInnerProductFactory
        = InnerProductFactory<TSpaces, InnerProductTerms..., NewTerm>;
      return NewInnerProductFactory{
               testSpaces,
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
      auto cFunc = Functions::detail::toGridViewFunction<GridView>(c);
      auto lhsBetaFunc
          = Functions::detail::toGridViewFunction<GridView>(lhsBeta);
      auto rhsBetaFunc
          = Functions::detail::toGridViewFunction<GridView>(rhsBeta);
      auto newTerm = make_IntegralTerm<lhsSpaceIndex,
                                       rhsSpaceIndex,
                                       integrationType,
                                       domainOfIntegration>
                                      (std::move(cFunc),
                                       std::move(lhsBetaFunc),
                                       std::move(rhsBetaFunc));
      using NewTerm = decltype(newTerm);
      using NewInnerProductFactory
        = InnerProductFactory<TSpaces, InnerProductTerms..., NewTerm>;
      return NewInnerProductFactory{
               testSpaces,
               std::tuple_cat(terms, std::make_tuple(std::move(newTerm)))
             };
    }

    InnerProduct<TSpaces, std::tuple<InnerProductTerms...>>
    create()
    {
      return {testSpaces, terms};
    }
  };

  template<class TSpaces>
  InnerProductFactory<TSpaces>
  innerProductWithSpace(const TSpaces& testSpaces)
  {
    return {testSpaces, std::tuple<>{}};
  }
}

#endif
