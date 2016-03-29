// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_WAVELET_SCATTERING_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_WAVELET_SCATTERING_HH

#include <tuple>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/dpg/assemble_types.hh>
#include <dune/dpg/assemble_helper.hh>
#include <dune/dpg/quadrature.hh>

#include <dune/istl/bvector.hh>

#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/functional/generation/make_fused_procedure.hpp>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Core>
#include <Eigen/SVD>

namespace Dune {


/**
 * \brief This constructs the right hand side vector of a DPG system.
 *
 * \tparam TestSpaces      tuple of test spaces
 * \tparam SolutionSpaces  tuple of solution spaces
 * \tparam FormulationType either SaddlepointFormulation or DPGFormulation
 */
template<class TestSpaces,
         class SolutionSpaces,
         class FormulationType>
class WaveletScatteringAssembler
{
private:
  class SVD {
    public:
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

      void applyToVector(Eigen::VectorXd& v) {
        v = kernelSVD.matrixU().leftCols(rank)
          * kernelSVD.singularValues().head(rank).asDiagonal()
          * kernelSVD.matrixV().leftCols(rank).adjoint() * v;
      }

      void setAccuracy(double accuracy) {
        using namespace Eigen;
        VectorXd singularValues = kernelSVD.singularValues();
        rank = singularValues.size();
        double err = 0;
        accuracy = accuracy * accuracy;
        while (err < accuracy && rank > 0) {
          rank -= 1;
          err += singularValues(rank) * singularValues(rank);
        }
        rank += 1;
        // TODO: If accuracy is low enough to allow rank = 0,
        //       this gives rank = 1.
      }

    private:
      Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::NoQRPreconditioner> kernelSVD;
      size_t rank;
  };

public:
  enum : unsigned int { dim = 2 };
  using Direction = FieldVector<double, dim>;

  WaveletScatteringAssembler () = delete;
  /**
   * \brief constructor for WaveletScatteringAssembler
   *
   * \note For your convenience, use make_WaveletScatteringAssembler() instead.
   */
  template<class Function>
  WaveletScatteringAssembler (const TestSpaces& testSpaces,
                              const SolutionSpaces& solutionSpaces,
                              const Function& kernel,
                              size_t num_s)
             : testSpaces(testSpaces),
               solutionSpaces(solutionSpaces),
               kernelSVD(kernel, num_s)
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
   * \param[in]  si          index of the scattering direction
   * \param[in]  accuracy    accuracy of the kernel approximation
   */
  template<size_t solutionSpaceIndex>
  void assembleScattering
         (BlockVector<FieldVector<double,1>>& scattering,
          const std::vector<BlockVector<FieldVector<double,1>>>& x,
          size_t si,
          double accuracy);

private:
  TestSpaces     testSpaces;
  SolutionSpaces solutionSpaces;
  SVD kernelSVD;
};

/**
 * \brief Creates a ScatteringAssembler for a DPG discretization,
 *        deducing the target type from the types of arguments.
 *
 * \param  testSpaces       a tuple of test spaces
 * \param  solutionSpaces   a tuple of solution spaces
 * \param  kernel           a function over (double, double) describing
 *                            the scattering kernel
 * \param  num_s            number of scattering directions
 */
template<class TestSpaces,
         class SolutionSpaces,
         class Function>
auto make_DPG_WaveletScatteringAssembler(
      const TestSpaces& testSpaces,
      const SolutionSpaces& solutionSpaces,
      const Function& kernel,
      size_t num_s)
     -> WaveletScatteringAssembler<TestSpaces, SolutionSpaces, DPGFormulation>
{
  return WaveletScatteringAssembler<TestSpaces, SolutionSpaces, DPGFormulation>
                            (testSpaces, solutionSpaces, kernel, num_s);
}


template<class TestSpaces,
         class SolutionSpaces,
         class FormulationType>
template<size_t solutionSpaceIndex>
void WaveletScatteringAssembler<TestSpaces, SolutionSpaces, FormulationType>::
assembleScattering(BlockVector<FieldVector<double,1> >& scattering,
                   const std::vector<BlockVector<FieldVector<double,1> >>& x,
                   size_t si, double accuracy)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  kernelSVD.setAccuracy(accuracy);

  constexpr bool isSaddlepoint =
        std::is_same<
             typename std::decay<FormulationType>::type
           , SaddlepointFormulation
        >::value;

  const int numS = x.size();

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
  auto testLocalView     = as_vector(transform(testSpaces, getLocalView()));
  auto solutionLocalView = as_vector(transform(solutionSpaces, getLocalView()));

  auto testLocalIndexSet     = as_vector(transform(testSpaces,
                                                   getLocalIndexSet()));
  auto solutionLocalIndexSet = as_vector(transform(solutionSpaces,
                                                   getLocalIndexSet()));

  for(const auto& e : elements(gridView)) {

    // Bind the local FE basis view to the current element
    /* TODO: only bind the space we use later */
    for_each(solutionLocalView, applyBind<decltype(e)>(e));
    for_each(testLocalView, applyBind<decltype(e)>(e));

    for_each(zip(solutionLocalIndexSet, solutionLocalView),
             make_fused_procedure(bindLocalIndexSet()));
    for_each(zip(testLocalIndexSet, testLocalView),
             make_fused_procedure(bindLocalIndexSet()));

    // Now get the local contribution to the right-hand side vector

    // Get the grid element from the local FE basis view
    typedef typename std::remove_reference<decltype(e)>::type Element;

    const int dim = Element::dimension;

    // Get set of shape functions for this element
    const auto& localFiniteElementTest =
        at_c<0>(testLocalView).tree().finiteElement();
    const auto& localFiniteElementSolution =
        at_c<solutionSpaceIndex>(solutionLocalView).tree().finiteElement();

    BlockVector<FieldVector<double,1> >
        localScattering(localFiniteElementTest.localBasis().size());
    localScattering = 0;

    /*TODO:
    - Find out what quadrature rule we are exactly using
    - Adapt quadrature to the kernel k also
    */
    int quadratureOrder = localFiniteElementSolution.localBasis().order()
                        + localFiniteElementTest.localBasis().order();
    /* Remark:
    - localFiniteElementTest.localBasis().order() is the degree of the polynomial
    - quad.size() is the number of quadrature points
    */

    const QuadratureRule<double, dim>& quad =
        QuadratureRules<double, dim>::rule(e.type(), quadratureOrder);

    // Loop over all quadrature points
    for ( size_t pt=0; pt < quad.size(); pt++ ) {

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();

      // The multiplicative factor in the integral transformation formula
      const double integrationElement = e.geometry().integrationElement(quadPos);

      // Evaluate all test shape function values at this quadrature point
      std::vector<FieldVector<double,1>> testShapeFunctionValues;
      localFiniteElementTest.localBasis()
          .evaluateFunction(quadPos, testShapeFunctionValues);

      Eigen::VectorXd uValues(numS);
      {
        std::vector<FieldVector<double,1>> shapeFunctionValues;
        localFiniteElementSolution.localBasis().
            evaluateFunction(quadPos, shapeFunctionValues);
        for( size_t scatteringAngle=0;
             scatteringAngle<numS; ++scatteringAngle) {
          double uValue = 0; // in direction of scatteringAngle
          // Evaluate all shape function values at this point
          std::vector<FieldVector<double,1>> shapeFunctionValues;
          localFiniteElementSolution.localBasis().
              evaluateFunction(quadPos, shapeFunctionValues);
          for (size_t j=0; j<shapeFunctionValues.size(); j++)
          {
            /* This assumes that solutionLocalIndexSet and
             * globalSolutionSpaceOffset don't change for different
             * scattering angles.
             */
            auto row =
                at_c<solutionSpaceIndex>(solutionLocalIndexSet).index(j)[0]
              + globalSolutionSpaceOffset;
            uValue += x[scatteringAngle][row] * shapeFunctionValues[j];
          }
          uValues(scatteringAngle) = uValue;
        }
      }
      kernelSVD.applyToVector(uValues);

      const double factor = uValues(si) * quad[pt].weight()
                            * integrationElement;
      for (size_t i=0, i_max=localScattering.size(); i<i_max; i++)
        localScattering[i] += factor * testShapeFunctionValues[i];
    }

    auto rhsCopier
        = localToGlobalRHSCopier
            <typename std::remove_reference<decltype(scattering)>::type>
            (scattering);
    rhsCopier(localScattering,
              at_c<0>(testLocalIndexSet),
              globalTestSpaceOffsets[0]
             );
  }
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_WAVELET_SCATTERING_HH
