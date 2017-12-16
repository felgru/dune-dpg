// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_QUADRATURE_HH
#define DUNE_DPG_QUADRATURE_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/quadraturerules/subsampledquadraturerule.hh>

#include "type_traits.hh"

namespace Dune {

namespace detail {
  template<class QuadratureRule>
  struct ChooseQuadratureImpl;

  template<int s, int dim>
  struct ChooseQuadratureImpl<Dune::SubsampledQuadratureRule<double, s, dim>>
  {
    template<class Element>
    static Dune::SubsampledQuadratureRule<double, s, dim>
    Quadrature(const Element& element,
               unsigned int quadratureOrder)
    {
      const Dune::QuadratureRule<double, dim>& quadSection =
            Dune::QuadratureRules<double, dim>::rule(element.type(),
                                                     quadratureOrder);
      SubsampledQuadratureRule<double, s, dim> quad(quadSection);
      return quad;
    }
  };

  template<int dim>
  struct ChooseQuadratureImpl<const Dune::QuadratureRule<double, dim>&>
  {
    template<class Element>
    static const Dune::QuadratureRule<double, dim>&
    Quadrature(const Element& element,
               unsigned int quadratureOrder)
    {
      const Dune::QuadratureRule<double, dim>& quad =
            Dune::QuadratureRules<double, dim>::rule(element.type(),
                                                     quadratureOrder);
      return quad;
    }
  };

  template<class LhsSpace, class RhsSpace, typename Element>
  struct ChooseQuadrature {
    static const int dim = Element::mydimension;
    static const bool useSubsampledQuadrature =
      is_SubsampledFiniteElement<LhsSpace>::value ||
      is_SubsampledFiniteElement<RhsSpace>::value;

    using Slhs = numberOfSamples<LhsSpace>;
    using Srhs = numberOfSamples<RhsSpace>;
    static constexpr int s =
        std::conditional<Slhs::value < Srhs::value, Srhs, Slhs>::type::value;

    using type
        = typename std::conditional
          < useSubsampledQuadrature
          , SubsampledQuadratureRule<double, s, dim>
          , const QuadratureRule<double, dim>&
          >::type;

    static type Quadrature(const Element& element,
                           unsigned int quadratureOrder)
    {
      return ChooseQuadratureImpl<type>
        ::template Quadrature(element, quadratureOrder);
    }
  };

} // End namespace detail

} // End namespace Dune

#endif
