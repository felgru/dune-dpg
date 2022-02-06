// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_SUBSAMPLED_QUADRATURE_RULE_HH
#define DUNE_GEOMETRY_SUBSAMPLED_QUADRATURE_RULE_HH

/** \file
 * \brief Construct summed quadrature rules from other quadrature rules
 */

#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dune/geometry/quadraturerules.hh>

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
     * \param quad Base quadrature rule.  Element type of this
     * rule must be cubical or triangle
     */
    SubsampledQuadratureRule(const Dune::QuadratureRule<ctype,dim>& quad)
      : QuadratureRule<ctype,dim>(quad.type(), quad.order())
    {
      if(quad.type().isCube()) {
        constexpr unsigned int numSections = power(s, dim);
        constexpr ctype volumeFraction = 1./static_cast<ctype>(numSections);

        this->reserve(numSections*quad.size());

        for (unsigned int i=0; i<numSections; i++) {
          auto alpha = multiindex(i);

          for (size_t q=0, qsize=quad.size(); q<qsize; q++) {

            Dune::FieldVector<ctype,dim> position = quad[q].position();
            for(unsigned int d=0; d<dim; d++) {
              position[d] += alpha[d];
              position[d] /= static_cast<ctype>(s);
            }
            this->emplace_back(position,
                               volumeFraction*quad[q].weight());

          }

        }

      } else if(quad.type().isTriangle()) {
        constexpr unsigned int numSections = power(s, dim);
        constexpr ctype volumeFraction = 1./static_cast<ctype>(numSections);

        this->reserve(numSections*quad.size());

        for (unsigned int line=0; line<s; ++line) {
          for(unsigned int i=0; i<2*(s-line)-1; ++i) {
            if(i%2) { // mirrored triangle
              for (size_t q=0, qsize=quad.size(); q<qsize; q++) {

                // mirror position accordingly
                Dune::FieldVector<ctype,dim> position = quad[q].position();
                ctype tmp = position[0];
                position[0] = 1 - position[1] + i/2;
                position[0] /= static_cast<ctype>(s);
                position[1] = 1 - tmp + line;
                position[1] /= static_cast<ctype>(s);

                this->emplace_back(position,
                                   volumeFraction*quad[q].weight());

              }
            } else {
              for (size_t q=0, qsize=quad.size(); q<qsize; q++) {

                Dune::FieldVector<ctype,dim> position = quad[q].position();
                position[0] += i/2;
                position[0] /= static_cast<ctype>(s);
                position[1] += line;
                position[1] /= static_cast<ctype>(s);

                this->emplace_back(position,
                                   volumeFraction*quad[q].weight());

              }
            }
          }
        }
      } else {
        DUNE_THROW(Dune::NotImplemented,
                   "SubsampledQuadratureRule only implemented on "
                   "cubes and triangles!");
      }

    }

  };

}

#endif   // DUNE_GEOMETRY_SUBSAMPLED_QUADRATURE_RULE_HH
