// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_FUNCTIONS_DISCRETEGLOBALBASISFUNCTION_HH
#define DUNE_DPG_FUNCTIONS_DISCRETEGLOBALBASISFUNCTION_HH

#include<type_traits>
#include <dune/dpg/functions/concepts.hh>
#include <dune/functions/functionspacebases/concepts.hh>

namespace Dune {

  namespace Functions {
    template<typename R, typename B, typename V>
    auto
    makeDiscreteGlobalBasisFunction(B&& basis, V&& vector);

    template<typename R, typename B, typename V>
    auto
    makeConstrainedDiscreteGlobalBasisFunction(const B& basis, const V& vector);
  }

  template<class FEBasis, class Vector,
      typename std::enable_if<models<Functions::Concept
                              ::GlobalBasis<typename FEBasis::GridView>,
                            FEBasis>()>::type* = nullptr>
  inline auto
  discreteGlobalBasisFunction(const FEBasis& feBasis, const Vector& u) {
    auto uFunction
        = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
              (feBasis, u);
    return uFunction;
  }

  template<class FEBasis, class Vector,
      typename std::enable_if<models<Functions::Concept::
            ConstrainedGlobalBasis<typename FEBasis::GridView>,
          FEBasis>()>::type* = nullptr>
  inline auto
  discreteGlobalBasisFunction(const FEBasis& feBasis, const Vector& u) {
    auto uFunction = Dune::Functions
        ::makeConstrainedDiscreteGlobalBasisFunction<double>(feBasis, u);
    return uFunction;
  }
} // end namespace Dune

#endif // DUNE_DPG_FUNCTIONS_DISCRETEGLOBALBASISFUNCTION_HH
