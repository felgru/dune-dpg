// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_APPROXIMATE_SCATTERING_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_APPROXIMATE_SCATTERING_HH

#include <memory>
#include <tuple>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/localevaluation.hh>
#include <dune/dpg/quadrature.hh>

#include <dune/istl/bvector.hh>

#include "svdkernelapproximation.hh"
#include "waveletkernelapproximation.hh"

namespace Dune {

/**
 * \brief This constructs the right hand side vector of a DPG system.
 *
 * \tparam SolutionSpace  the GlobalBasis of the solution spaces (on host grid)
 */
template<class SolutionSpace,
         class KernelApproximation>
class ApproximateScatteringAssembler
{
public:
  enum : unsigned int { dim = SolutionSpace::GridView::dimension };
  using Direction = FieldVector<double, dim>;

  using GridData = BlockVector<FieldVector<double,1>>;

  ApproximateScatteringAssembler () = delete;
  /**
   * \brief constructor for ApproximateScatteringAssembler
   *
   * \param solutionSpace    a reference to a solution space (on the host grid)
   * \param kernel           an approximation of the scattering kernel
   *
   * \note For your convenience, use make_ApproximateScatteringAssembler()
   *       instead.
   */
  ApproximateScatteringAssembler (
      const SolutionSpace& solutionSpace,
      const KernelApproximation& kernel)
             : solutionSpace(solutionSpace),
               kernelApproximation(kernel)
  {}

  /**
   * \brief Assemble a vector representing the scattering integral
   *        in the solution space basis for a given set of
   *        discrete ordinate solutions.
   *
   * \todo For now, the scattering kernel is assumed to be constant
   *       in space and normed to integral 1.
   *
   * \param[out] scattering  the scattering vectors
   * \param[in]  x           the vectors of the solutions of the
   *                           previous iteration
   */
  void computeScattering
          (std::vector<GridData>& scattering,
           const std::vector<GridData>& x) const;

private:
  const SolutionSpace& solutionSpace;
  const KernelApproximation& kernelApproximation;
};

/**
 * \brief Creates a ScatteringAssembler for a DPG discretization,
 *        deducing the target type from the types of arguments.
 *
 * \param  solutionSpace    a reference to a solution spaces (on host grid).
 *                          It is assumed that the solutions given later have
 *                          all been interpolated to this solutionSpace.
 * \param  kernel           an approximation of the scattering kernel
 */
template<class SolutionSpace,
         class KernelApproximation>
auto make_ApproximateScatteringAssembler(
      const SolutionSpace& solutionSpace,
      const KernelApproximation& kernel)
     -> ApproximateScatteringAssembler<SolutionSpace, KernelApproximation>
{
  return ApproximateScatteringAssembler<SolutionSpace, KernelApproximation>
                            (solutionSpace, kernel);
}


template<class SolutionSpace,
         class KernelApproximation>
void ApproximateScatteringAssembler<SolutionSpace, KernelApproximation>::
computeScattering(std::vector<GridData>& scattering,
                     const std::vector<GridData>& x) const
{
  using namespace Dune::detail;

  const size_t numDoFs = x[0].size();
  const size_t numSscattered = scattering.size();
  for(size_t i = 0; i < numSscattered; i++) {
    scattering[i].resize(numDoFs);
  }

  const size_t numS = x.size();
  for(size_t i = 0; i < numS; i++) {
    assert(x[i].size() == numDoFs);
  }
  for(size_t row = 0; row < numDoFs; row++) {
    Eigen::VectorXd uValues(numS);

    // get values for all directions at position x_row
    for (size_t scatteringAngle=0;
         scatteringAngle<numS; ++scatteringAngle) {
      uValues(scatteringAngle) = x[scatteringAngle][row];
    }

    kernelApproximation.applyToVector(uValues);

    for (size_t scatteringAngle=0;
         scatteringAngle<numSscattered; ++scatteringAngle) {
      scattering[scatteringAngle][row] = uValues(scatteringAngle);
    }
  }
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_APPROXIMATE_SCATTERING_HH
