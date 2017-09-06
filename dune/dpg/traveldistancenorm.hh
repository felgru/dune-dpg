// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_TRAVELDISTANCENORM_HH
#define DUNE_DPG_TRAVELDISTANCENORM_HH

#include <cfloat>
#include <cmath>

namespace Dune {

namespace detail {
  // This function expects data transformed to the reference triangle:
  // A point in the inflow boundary of the reference cell and the
  // transport direction transformed under the global-to-local mapping.
  template<class ReferenceCellCoordinate, class ReferenceCellDirection>
  double travelDistance(
      const ReferenceCellCoordinate& x,
      const ReferenceCellDirection& beta)
  {
    if(x[0]) { /* x[0] != 0 */
      if(x[1]) { /* x[1] != 0 */
        const double a = (beta[0]!=0.)?(x[1] - x[0]*beta[1]/beta[0])
                                      :(-DBL_MAX);
        if(0 <= a && a <= 1)
          return -x[0]/beta[0];
        else
          return -x[1]/beta[1];
      } else { /* x[1] == 0 */
        const double a = (beta[0]!=0.)?(-beta[1]/beta[0]*x[0]):(-DBL_MAX);
        if(0 <= a && a <= 1)
          return -x[0]/beta[0];
        else
          return (1-x[0])/(beta[0]+beta[1]);
      }
    } else { /* x[0] == 0 */
      const double b = (beta[1]!=0.)?(-beta[0]/beta[1]*x[1]):(-DBL_MAX);
      if(0 <= b && b <= 1)
        return -x[1]/beta[1];
      else
        return (1-x[1])/(beta[0]+beta[1]);
    }
  }

  template<class FaceGeometryInElement, class ReferenceCellDirection>
  double splitPointOfInflowFaceInTriangle(
      const FaceGeometryInElement& faceGeometryInElement,
      ReferenceCellDirection& referenceBeta)
  {
    // This gets a bit ugly as we have to check the orientation of the face
    double splitPoint;
    const double tol = 1e-8;
    // Take the 1-norm as it is cheaper to compute.
    const double referenceBetaNorm = fabs(referenceBeta[0])
                                   + fabs(referenceBeta[1]);
    // We use the relative error w.r.t. the norm of referenceBeta as
    // the numerical errors in the computation of referenceBeta seem
    // to scale with decreasing diameter of the element which is
    // proportional to an increase in the norm of referenceBeta.
    // When using an absolute error instead, we observed that tol was
    // always too small for later iterations in the grid refinement.
    if(fabs(referenceBeta[0]) < tol * referenceBetaNorm) referenceBeta[0] = 0.;
    if(fabs(referenceBeta[1]) < tol * referenceBetaNorm) referenceBeta[1] = 0.;

    auto corner = faceGeometryInElement.global({0});
    if(referenceBeta[0] > 0) {
      if((corner
          - FieldVector<double,2>{0.,0.}).two_norm() < tol)
        splitPoint = -referenceBeta[1]/referenceBeta[0];
      else
        splitPoint = 1.+referenceBeta[1]/referenceBeta[0];
    } else if(referenceBeta[1] > 0) {
      if((corner
          - FieldVector<double,2>{0.,0.}).two_norm() < tol)
        splitPoint = -referenceBeta[0]/referenceBeta[1];
      else
        splitPoint = 1.+referenceBeta[0]/referenceBeta[1];
    } else {
      if((corner
          - FieldVector<double,2>{0.,1.}).two_norm() < tol)
        splitPoint = referenceBeta[0]/(referenceBeta[0]+referenceBeta[1]);
      else
        splitPoint = 1.-referenceBeta[0]/(referenceBeta[0]+referenceBeta[1]);
    }
    if(fabs(splitPoint)< tol)         splitPoint = 0.;
    else if(fabs(splitPoint-1.)< tol) splitPoint = 1.;

    assert(splitPoint >= 0 && splitPoint <= 1);
    return splitPoint;
  }

  template <class Geometry, class SubGeometry, int dim>
  FieldVector<double,dim> referenceBeta(
      const Geometry& geometry,
      const SubGeometry& subGeometryInReferenceElement,
      const FieldVector<double, dim>& beta)
  {
    static_assert(dim==2, "Computation of transport direction on reference"
                          " cell only implemented in 2d!");
    /* This won't work for curvilinear elements, but they don't seem
     * to be supported by UG anyway. */
    const auto& jacobianSubInverse
        = subGeometryInReferenceElement.jacobianInverseTransposed({0., 0.});
    const auto& jacobianInverse = geometry.jacobianInverseTransposed({0., 0.});
    FieldVector<double,dim> referenceBetaSub, referenceBeta;
    jacobianInverse.mtv(beta, referenceBeta);
    jacobianSubInverse.mtv(referenceBeta, referenceBetaSub);

    return referenceBetaSub;
  }
} // end namespace detail

} // end namespace Dune

#endif
