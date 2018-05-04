// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_LINEARFORMFACTORY_HH
#define DUNE_DPG_LINEARFORMFACTORY_HH

#include <tuple>
#include <type_traits>

#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/functions/gridviewfunctions.hh>
#include <dune/dpg/linearform.hh>
#include <dune/dpg/linearfunctionalterm.hh>
#include <dune/dpg/linearintegralterm.hh>
#include <dune/dpg/spacetuple.hh>

namespace Dune {
  template<class TSpaces, class... LinearFormTerms>
  class LinearFormFactory;

  template<class TSpaces>
  LinearFormFactory<TSpaces>
  linearFormWithSpace(const TSpaces& testSpaces);

  template<class TSpaces, class... LinearFormTerms>
  class LinearFormFactory {
    static_assert(is_SpaceTuplePtr<TSpaces>::value,
        "TSpaces needs to be a SpaceTuplePtr!");
    using TestSpacesPtr = TSpaces;
    using TestSpaces = typename TestSpacesPtr::element_type;
    using LinearFormTermsTuple = std::tuple<LinearFormTerms...>;
    using GridView = typename std::tuple_element_t<0,TestSpaces>::GridView;

    LinearFormFactory(const TestSpacesPtr& testSpaces,
                        const LinearFormTermsTuple& terms)
      : testSpaces(testSpaces)
      , terms(terms)
      {}

    friend
    LinearFormFactory<TSpaces>
    linearFormWithSpace<TSpaces>(const TSpaces& testSpaces);

    template<class TSpaces_, class... LinearFormTerms_>
    friend class LinearFormFactory;

    TestSpacesPtr testSpaces;
    LinearFormTermsTuple terms;

  public:
    template<size_t spaceIndex,
             LinearIntegrationType integrationType,
             DomainOfIntegration domainOfIntegration,
             class Factor>
    auto addIntegralTerm(Factor c)
    {
      auto cFunc = Functions::detail::toGridViewFunction<GridView>(c);
      auto newTerm = make_LinearIntegralTerm<spaceIndex,
                                             integrationType,
                                             domainOfIntegration>
                                                    (std::move(cFunc));
      using NewTerm = decltype(newTerm);
      using NewLinearFormFactory
        = LinearFormFactory<TSpaces, LinearFormTerms..., NewTerm>;
      return NewLinearFormFactory{
               testSpaces,
               std::tuple_cat(terms, std::make_tuple(std::move(newTerm)))
             };
    }

    template<size_t spaceIndex,
             LinearIntegrationType integrationType,
             DomainOfIntegration domainOfIntegration,
             class Factor,
             class Direction>
    auto addIntegralTerm(Factor c, Direction beta)
    {
      auto cFunc = Functions::detail::toGridViewFunction<GridView>(c);
      auto betaFunc = Functions::detail::toGridViewFunction<GridView>(beta);
      auto newTerm = make_LinearIntegralTerm<spaceIndex,
                                             integrationType,
                                             domainOfIntegration>
                              (std::move(cFunc), std::move(betaFunc));
      using NewTerm = decltype(newTerm);
      using NewLinearFormFactory
        = LinearFormFactory<TSpaces, LinearFormTerms..., NewTerm>;
      return NewLinearFormFactory{
               testSpaces,
               std::tuple_cat(terms, std::make_tuple(std::move(newTerm)))
             };
    }

    template<size_t spaceIndex,
             class SolutionSpace,
             class FunctionalVector>
    auto addFunctionalTerm(const FunctionalVector& functionalVector,
                           const SolutionSpace& solutionSpace)
    {
      auto newTerm = make_LinearFunctionalTerm<spaceIndex, SolutionSpace,
                                               FunctionalVector>
                                              (functionalVector, solutionSpace);
      using NewTerm = decltype(newTerm);
      using NewLinearFormFactory
        = LinearFormFactory<TSpaces, LinearFormTerms..., NewTerm>;
      return NewLinearFormFactory{
               testSpaces,
               std::tuple_cat(terms, std::make_tuple(std::move(newTerm)))
             };
    }

    template<size_t spaceIndex,
             IntegrationType integrationType,
             class SolutionSpace,
             class FunctionalVector,
             class Factor, class Direction,
             typename std::enable_if<
                         integrationType == IntegrationType::normalVector
                      || integrationType ==
                                  IntegrationType::travelDistanceWeighted
                      >::type*
               = nullptr
            >
    auto addSkeletalFunctionalTerm(
        const FunctionalVector& functionalVector,
        const SolutionSpace& solutionSpace,
        Factor c, Direction beta)
    {
      auto cFunc = Functions::detail::toGridViewFunction<GridView>(c);
      auto betaFunc = Functions::detail::toGridViewFunction<GridView>(beta);
      auto newTerm = make_SkeletalLinearFunctionalTerm
                      <spaceIndex, integrationType>
                      (functionalVector, solutionSpace,
                       std::move(cFunc), std::move(betaFunc));
      using NewTerm = decltype(newTerm);
      using NewLinearFormFactory
        = LinearFormFactory<TSpaces, LinearFormTerms..., NewTerm>;
      return NewLinearFormFactory{
               testSpaces,
               std::tuple_cat(terms, std::make_tuple(std::move(newTerm)))
             };
    }

    LinearForm<TSpaces, std::tuple<LinearFormTerms...>>
    create()
    {
      return {testSpaces, terms};
    }
  };

  template<class TSpaces>
  LinearFormFactory<TSpaces>
  linearFormWithSpace(const TSpaces& testSpaces)
  {
    return {testSpaces, std::tuple<>{}};
  }
}

#endif
