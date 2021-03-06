// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_SPLIT_QUADRATURE_RULE_HH
#define DUNE_GEOMETRY_SPLIT_QUADRATURE_RULE_HH

/** \file
 * \brief Construct summed quadrature rules from other quadrature rules
 */

#include <dune/common/fvector.hh>

#include <dune/geometry/quadraturerules.hh>

namespace Dune {

  /** \brief Construct summed quadrature rules from other quadrature rules
   *
   * \tparam ctype Type used for coordinates and quadrature weights
   */
  template <class ctype>
  class SplitQuadratureRule
      : public Dune::QuadratureRule<ctype,1>
  {
    public:
    /** \brief Construct summed quadrature rule
     * \param quad Base quadrature rule.
     * \param splitPoint point between 0 and 1 where the quadrature is split
     */
    SplitQuadratureRule(const Dune::QuadratureRule<ctype,1>& quad,
        ctype splitPoint)
      : QuadratureRule<ctype,1>(quad.type(), quad.order())
    {
      this->reserve(2*quad.size());
      ctype volumeFraction = splitPoint;
      for (const auto& q : quad) {
        const Dune::FieldVector<ctype,1> position
            = {q.position()[0] * volumeFraction + 0.};
        this->emplace_back(position,
                           volumeFraction * q.weight());
      }
      volumeFraction = 1 - splitPoint;
      for (const auto& q : quad) {
        const Dune::FieldVector<ctype,1> position
            = {q.position()[0] * volumeFraction + splitPoint};
        this->emplace_back(position,
                           volumeFraction * q.weight());

      }
    }

  };

}

#endif   // DUNE_GEOMETRY_SPLIT_QUADRATURE_RULE_HH
