// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_SUBSAMPLED_QUADRATURE_RULE_HH
#define DUNE_GEOMETRY_SUBSAMPLED_QUADRATURE_RULE_HH

/** \file
 * \brief Construct summed quadrature rules from other quadrature rules
 */

#include <dune/common/fvector.hh>
#include <dune/common/power.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/virtualrefinement.hh>

namespace Dune {

  /** \brief Construct summed quadrature rules from other quadrature rules
   *
   * \tparam ctype Type used for coordinates and quadrature weights
   * \tparam s number of subsamples
   * \tparam dim Dimension of the reference element
   */
  template <class ctype, int s, int dim>
  class SubsampledQuadratureRule
      : public Dune::QuadratureRule<ctype,dim>
  {
    private:
    // Return i as a d-digit number in the s-ary system
    static Dune::FieldVector<int,dim> multiindex (int i)
    {
      Dune::FieldVector<int,dim> alpha;
      for (int j=0; j<dim; j++)
      {
        alpha[j] = i % s;
        i = i/s;
      }
      return alpha;
    }

    public:
    /** \brief Construct summed quadrature rule
     * \param quad Base quadrature rule.  Element type of this rule must be simplex
     */
    SubsampledQuadratureRule(const Dune::QuadratureRule<ctype,dim>& quad)
      : QuadratureRule<ctype,dim>(quad.type(), quad.order())
    {
      // Currently only works for cubes
      assert(quad.type().isCube());

      unsigned int numSections = StaticPower<s,dim>::power;
      ctype volumeFraction = 1./(ctype)numSections;

      this->reserve(numSections*quad.size());

      for (unsigned int i=0; i<numSections; i++) {
        auto alpha = multiindex(i);

        for (size_t q=0; q<quad.size(); q++) {

          Dune::FieldVector<ctype,dim> position = quad[q].position();
          for(unsigned int d=0; d<dim; d++) {
            position[d] += alpha[d];
            position[d] /= (ctype)s;
          }
          // TODO: use emplace_back
          this->push_back(Dune::QuadraturePoint<ctype,dim>(position,
                      volumeFraction*quad[q].weight()));

        }

      }

    }

  };

}

#endif   // DUNE_GEOMETRY_SUBSAMPLED_QUADRATURE_RULE_HH
