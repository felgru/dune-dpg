// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_APPROXIMATE_SCATTERING_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_APPROXIMATE_SCATTERING_HH

#include <tuple>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/localevaluation.hh>
#include <dune/dpg/quadrature.hh>

#include <dune/istl/bvector.hh>

#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/functional/generation/make_fused_procedure.hpp>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Core>
#include <Eigen/SVD>

namespace Dune {

namespace ScatteringKernelApproximation {
  class SVD {
    public:
      enum : unsigned int { dim = 2 };
      using Direction = FieldVector<double, dim>;

      SVD() = delete;
      SVD(const SVD&) = delete;

      template<class Function>
      SVD(const Function& kernel, size_t num_s)
        : kernelSVD(num_s, num_s, Eigen::ComputeThinU | Eigen::ComputeThinV),
          rank(num_s) {
        using namespace Eigen;
        using namespace boost::math::constants;
        MatrixXd kernelMatrix(num_s, num_s);
        for(size_t j = 0; j < num_s; ++j) {
          Direction s_j = {cos(2*pi<double>()*j/num_s),
                           sin(2*pi<double>()*j/num_s)};
          for(size_t i = 0; i < num_s; ++i) {
            Direction s_i = {cos(2*pi<double>()*i/num_s),
                             sin(2*pi<double>()*i/num_s)};
            // TODO: maybe use a higher order quadrature
            kernelMatrix(i,j) = kernel(s_i, s_j)/(num_s*num_s);
          }
        }
        /* initialize SVD of kernel (using Eigen) */
        kernelSVD.compute(kernelMatrix);
      }

      void applyToVector(Eigen::VectorXd& v) const {
        v = kernelSVD.matrixU().leftCols(rank)
          * kernelSVD.singularValues().head(rank).asDiagonal()
          * kernelSVD.matrixV().leftCols(rank).adjoint() * v;
      }

      void setAccuracy(double accuracy) {
        using namespace Eigen;
        VectorXd singularValues = kernelSVD.singularValues();
        size_t i = singularValues.size() - 1;
        double err = 0,
               rank_err = singularValues(i) * singularValues(i);
        accuracy = accuracy * accuracy;
        while (err + rank_err < accuracy && i > 0) {
          err += rank_err;
          i -= 1;
          rank_err = singularValues(i) * singularValues(i);
        }
        rank = i+1;
        // TODO: If accuracy is low enough to allow rank = 0,
        //       this gives rank = 1.
      }

    private:
      Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::NoQRPreconditioner> kernelSVD;
      size_t rank;
  };
}

/**
 * \brief This constructs the right hand side vector of a DPG system.
 *
 * \tparam TestSpaces      tuple of test spaces
 * \tparam SolutionSpaces  tuple of solution spaces
 * \tparam FormulationType either SaddlepointFormulation or DPGFormulation
 */
template<class TestSpaces,
         class SolutionSpaces,
         class KernelApproximation,
         class FormulationType>
class ApproximateScatteringAssembler
{
public:
  enum : unsigned int { dim = 2 };
  using Direction = FieldVector<double, dim>;

  ApproximateScatteringAssembler () = delete;
  /**
   * \brief constructor for ApproximateScatteringAssembler
   *
   * \param si  index of the scattering direction
   *
   * \note For your convenience, use make_ApproximateScatteringAssembler()
   *       instead.
   */
  ApproximateScatteringAssembler (const TestSpaces& testSpaces,
                                  const SolutionSpaces& solutionSpaces,
                                  const KernelApproximation& kernel,
                                  size_t si)
             : testSpaces(testSpaces),
               solutionSpaces(solutionSpaces),
               si(si),
               kernelApproximation(kernel)
  {}

  /**
   * \brief Assemble the vector corresponding to the scattering integral
   * for a given set of discrete ordinate solutions.
   *
   * \todo: For now, the scattering kernel is assumed to be constant
   *        and normed to integral 1.
   *
   * \param[out] scattering  the scattering vector
   * \param[in]  x           the vectors of the solutions of the
   *                           previous iteration
   */
  template<size_t solutionSpaceIndex>
  void assembleScattering
         (BlockVector<FieldVector<double,1>>& scattering,
          const std::vector<BlockVector<FieldVector<double,1>>>& x);

private:
  TestSpaces     testSpaces;
  SolutionSpaces solutionSpaces;
  const size_t si;
  const KernelApproximation& kernelApproximation;
};

/**
 * \brief Creates a ScatteringAssembler for a DPG discretization,
 *        deducing the target type from the types of arguments.
 *
 * \param  testSpaces       a tuple of test spaces
 * \param  solutionSpaces   a tuple of solution spaces
 * \param  kernel           an approximation of the scattering kernel
 * \param  si               index of scattering direction (0 < si < numS)
 */
template<class TestSpaces,
         class SolutionSpaces,
         class KernelApproximation>
auto make_DPG_ApproximateScatteringAssembler(
      const TestSpaces& testSpaces,
      const SolutionSpaces& solutionSpaces,
      const KernelApproximation& kernel,
      size_t si)
     -> ApproximateScatteringAssembler<TestSpaces, SolutionSpaces, KernelApproximation, DPGFormulation>
{
  return ApproximateScatteringAssembler<TestSpaces, SolutionSpaces, KernelApproximation, DPGFormulation>
                            (testSpaces, solutionSpaces, kernel, si);
}

namespace detail {
  template <class TestSpace,
            class SolutionSpace,
            bool = is_RefinedFiniteElement<TestSpace>::value,
            bool = is_RefinedFiniteElement<SolutionSpace>::value>
  struct GetLocalApproximateScattering;
}
#include "approximatelocalscattering_uu_impl.hh"
#include "approximatelocalscattering_ru_impl.hh"


template<class TestSpaces,
         class SolutionSpaces,
         class KernelApproximation,
         class FormulationType>
template<size_t solutionSpaceIndex>
void ApproximateScatteringAssembler<TestSpaces, SolutionSpaces, KernelApproximation, FormulationType>::
assembleScattering(BlockVector<FieldVector<double,1> >& scattering,
                   const std::vector<BlockVector<FieldVector<double,1> >>& x)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  constexpr bool isSaddlepoint =
        std::is_same<
             typename std::decay<FormulationType>::type
           , SaddlepointFormulation
        >::value;

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  size_t globalSolutionSpaceOffsets[std::tuple_size<SolutionSpaces>::value];
  size_t globalTotalTestSize = 0;

  globalTotalTestSize =
    fold(zip(globalTestSpaceOffsets, testSpaces),
         (size_t)0, globalOffsetHelper());

  if(!isSaddlepoint)
  {
    for(size_t i=0; i<std::tuple_size<TestSpaces>::value; ++i)
    {
      globalTestSpaceOffsets[i] = 0;
    }
  }

  size_t globalTotalSolutionSize =
    fold(zip(globalSolutionSpaceOffsets, solutionSpaces),
         isSaddlepoint?globalTotalTestSize:0, globalOffsetHelper());
  globalTotalSolutionSize -= globalSolutionSpaceOffsets[0];

  double globalSolutionSpaceOffset =
      globalSolutionSpaceOffsets[solutionSpaceIndex];

  scattering.resize(globalTotalTestSize
                    + (isSaddlepoint?globalTotalSolutionSize:0));
  scattering = 0;

  // Views on the FE bases on a single element
  auto testLocalViews     = as_vector(transform(testSpaces, getLocalView()));
  auto solutionLocalViews = as_vector(transform(solutionSpaces, getLocalView()));

  auto testLocalIndexSets     = as_vector(transform(testSpaces,
                                                    getLocalIndexSet()));
  auto solutionLocalIndexSets = as_vector(transform(solutionSpaces,
                                                    getLocalIndexSet()));

  for(const auto& e : elements(gridView)) {

    // Bind the local FE basis view to the current element
    /* TODO: only bind the space we use later */
    for_each(solutionLocalViews, applyBind<decltype(e)>(e));
    for_each(testLocalViews, applyBind<decltype(e)>(e));

    for_each(zip(solutionLocalIndexSets, solutionLocalViews),
             make_fused_procedure(bindLocalIndexSet()));
    for_each(zip(testLocalIndexSets, testLocalViews),
             make_fused_procedure(bindLocalIndexSet()));

    // Now get the local contribution to the right-hand side vector

    BlockVector<FieldVector<double,1> >
        localScattering(at_c<0>(testLocalViews).size());
    localScattering = 0;

    detail::GetLocalApproximateScattering
      < std::tuple_element_t<0, TestSpaces>
      , std::tuple_element_t<solutionSpaceIndex, SolutionSpaces>
      >::interiorImpl(
            at_c<0>(testLocalViews),
            at_c<solutionSpaceIndex>(solutionLocalViews),
            localScattering,
            at_c<solutionSpaceIndex>(solutionLocalIndexSets),
            globalSolutionSpaceOffset,
            // unsigned int quadratureOrder,
            e,
            kernelApproximation,
            x,
            si);

    auto rhsCopier
        = localToGlobalRHSCopier
            <typename std::remove_reference<decltype(localScattering)>::type,
             typename std::remove_reference<decltype(scattering)>::type>
            (localScattering, scattering);
    // TODO: We should probably not hard-code the test space index
    //       and testLocalOffset.
    rhsCopier(at_c<0>(testLocalViews),
              at_c<0>(testLocalIndexSets),
              0,
              globalTestSpaceOffsets[0]
             );
  }
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_APPROXIMATE_SCATTERING_HH
