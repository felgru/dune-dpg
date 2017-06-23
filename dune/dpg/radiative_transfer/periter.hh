// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH

#include <chrono>
#include <set>
#include <vector>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/istl/umfpack.hh>

#include <dune/functions/functionspacebases/hangingnodep2nodalbasis.hh>
#include <dune/functions/functionspacebases/pqkdgrefineddgnodalbasis.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pqktracenodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/pqksubsampleddgbasis.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/dpg/boundarytools.hh>
#include <dune/dpg/errortools.hh>
#include <dune/dpg/functionplotter.hh>
#include <dune/dpg/functions/interpolate.hh>
#include <dune/dpg/functions/subgridinterpolation.hh>
#include <dune/dpg/functions/refinementinterpolation.hh>
#include <dune/dpg/linearfunctionalterm.hh>
#include <dune/dpg/radiative_transfer/approximate_scattering.hh>
#include <dune/dpg/radiative_transfer/subgridprojection.hh>
#include <dune/dpg/rhs_assembler.hh>
#include <dune/dpg/dpg_system_assembler.hh>
#include <dune/dpg/type_traits.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <dune/subgrid/subgrid.hh>
#pragma GCC diagnostic pop

#include <boost/math/constants/constants.hpp>

namespace Dune {

template<class SubGrid>
std::set<typename SubGrid::HostGridType::GlobalIdSet::IdType>
saveSubGridToIdSet(const SubGrid& subGrid)
{
  auto& idSet = subGrid.getHostGrid().globalIdSet();
  auto subGridView = subGrid.leafGridView();
  std::set<typename std::decay_t<decltype(idSet)>::IdType> subGridElements;
  for(const auto& e : elements(subGridView)) {
    subGridElements.insert(idSet.id(subGrid.template getHostEntity<0>(e)));
  }
  return subGridElements;
}

template<class SubGrid>
std::unique_ptr<SubGrid>
restoreSubGridFromIdSet(
    std::set<typename SubGrid::HostGridType::GlobalIdSet::IdType>&& idSet,
    typename SubGrid::HostGridType& hostGrid)
{
  auto subGrid = std::make_unique<SubGrid>(hostGrid);
  subGrid->createBegin();
  subGrid->insertSet(idSet);
  subGrid->createEnd();
  subGrid->setMaxLevelDifference(1);
  return subGrid;
}

template<class SubGrid>
std::unique_ptr<SubGrid>
restoreSubGridFromIdSet(
    std::set<typename SubGrid::HostGridType::GlobalIdSet::IdType>& idSet,
    typename SubGrid::HostGridType& hostGrid)
{
  std::set<typename SubGrid::HostGridType::GlobalIdSet::IdType>
    idSetCopy(idSet);
  return restoreSubGridFromIdSet<SubGrid>(std::move(idSetCopy), hostGrid);
}

enum class PlotSolutions {
  doNotPlot,
  plotOuterIterations,
  plotLastIteration
};

/**
 * This class describes the Periter algorithm for radiative transfer problems
 *
 * \tparam ScatteringKernelApproximation
 *         specifies the method used to approximate the scattering kernel
 * \tparam RHSApproximation  if right hand side and lifting of boundary
 *                           values are finite element functions, set this
 *                           to FeRHSandBoundary, otherwise set this to
 *                           ApproximateRHSandBoundary
 */
template<class ScatteringKernelApproximation, class RHSApproximation>
class Periter {
  public:

  /**
   * Solve a radiative transfer problem using the Periter algorithm
   *
   * \param hostGrid
   * \param f  right hand side function
   * \param g  lifting of the boundary values
   * \param gDeriv  derivative of g in direction s
   * \param sigma
   * \param kernel  the scattering kernel, e.g. a Henyey–Greenstein kernel
   * \param numS  number of directions used in the discretization
   * \param rho  the contraction parameter ρ
   * \param CT  the constant C_T from the paper
   * \param targetAccuracy  periter solves up to this accuracy
   * \param maxNumberOfIterations  ... or up to the given number of iterations
   *                               (whatever comes first)
   * \param plotSolutions  specifies when to create .vtu files for plotting
   *                       the solution
   */
  template<class Grid, class F, class G, class GDeriv, class Kernel>
  void solve(Grid& hostGrid,
             const F& f,
             const G& g,
             const GDeriv& gDeriv,
             double sigma,
             const Kernel& kernel,
             unsigned int numS,
             double rho,
             double CT,
             double targetAccuracy,
             unsigned int maxNumberOfIterations,
             PlotSolutions plotSolutions = PlotSolutions::doNotPlot);

  private:
  using VectorType = BlockVector<FieldVector<double,1>>;

  /**
   * Apply the scattering integral to a solution x
   *
   * This corresponds to [K, u_n, κ_1 * η] in the Periter algorithm
   * (see Dahmen, Gruber, Mula).
   *
   * \param kernelApproximation an approximation to the scattering kernel
   * \param x  solution to which we want to apply the scattering kernel
   * \param solutionSpaces
   * \param accuracy
   */
  template<class SolutionSpaces, class HostGrid>
  static std::vector<VectorType> apply_scattering(
      ScatteringKernelApproximation& kernelApproximation,
      const std::vector<VectorType>& x,
      const std::vector<std::shared_ptr<SolutionSpaces>>& solutionSpaces,
      const HostGrid& hostGrid,
      double accuracy);
};

struct FeRHSandBoundary {};
struct ApproximateRHSandBoundary {};

#ifndef DOXYGEN
namespace detail {
  constexpr size_t dim = 2;
  using VectorType = BlockVector<FieldVector<double,1>>;
  using Direction = FieldVector<double, dim>;

  template<class FEHostBasis, class Grids,
           class F, class G, class GDeriv>
  inline void approximate_rhs_and_bv (
      std::vector<VectorType>& rhsFunctional,
      Grids&,
      double,
      const FEHostBasis& hostGridBasis,
      const std::vector<Direction>& sVector,
      F& f, G& g, GDeriv& gDeriv, double sigma,
      FeRHSandBoundary) {
    static_assert(!is_RefinedFiniteElement<FEHostBasis>::value,
        "Functions::interpolate won't work for refined finite elements");
    size_t numS = sVector.size();
    for(unsigned int i = 0; i < numS; ++i)
    {
      VectorType gInterpolation(hostGridBasis.size());
      const Direction s = sVector[i];
      Functions::interpolate(hostGridBasis, gInterpolation,
          [s,&f,&g,&gDeriv,sigma](const Direction& x)
          { return f(x,s) + gDeriv(x,s) - sigma*g(x,s); });

      // Add gInterpolate to first hostGridBasis.size() entries of
      // rhsFunctional[i].
      using Iterator = std::decay_t<decltype(rhsFunctional[i].begin())>;
      for(Iterator rIt=rhsFunctional[i].begin(),
                   rEnd=rhsFunctional[i].begin()+hostGridBasis.size(),
                   gIt=gInterpolation.begin(); rIt!=rEnd; ++rIt, ++gIt) {
        *rIt += *gIt;
      }
    }
  }


  // Refine grid until accuracy kappa2*eta is reached for
  // approximation of rhs and boundary values.
  template<class FEBases, class Grids,
           class F, class G, class GDeriv>
  inline void approximate_rhs_and_bv (
      std::vector<VectorType>& rhsFunctional,
      Grids& grids,
      double accuracy,
      const std::vector<std::shared_ptr<FEBases>>& solutionSpaces,
      const std::vector<Direction>& sVector,
      F& f, G& g, GDeriv& gDeriv, double sigma,
      ApproximateRHSandBoundary) {
    DUNE_THROW(Dune::NotImplemented,
        "Implementation of approximate_rhs_and_bv for non-FE right hand side "
        "functions currently broken.");
    using Grid = typename Grids::value_type::element_type;
    using EntitySeed = typename Grid::template Codim<0>::Entity::EntitySeed;
    size_t numS = sVector.size();

    const size_t spaceIndex = 0;
    using FEBasisInterior = std::tuple_element_t<spaceIndex, FEBases>;
    static_assert(!is_RefinedFiniteElement<FEBasisInterior>::value,
        "Functions::interpolate won't work for refined finite elements");

    std::vector<VectorType> boundaryValues(numS);
    for(unsigned int i = 0; i < numS; ++i)
    {
      auto& feBasisInterior = std::get<spaceIndex>(*solutionSpaces[i]);
      const unsigned int maxNumberOfRefinements = 3;
      for(unsigned int refinement = 1; ; refinement++) {
        {
          VectorType gInterpolation(feBasisInterior.size());
          const Direction s = sVector[i];
          Functions::interpolate(feBasisInterior, gInterpolation,
              [s,&f,&g,&gDeriv,sigma](const Direction& x)
              { return f(x,s) + gDeriv(x,s) - sigma*g(x,s); });

          std::swap(boundaryValues[i], gInterpolation);
        }
        const Direction s = sVector[i];
        auto gExact = Functions::makeGridViewFunction(
              [s,&f,&g,&gDeriv,sigma](const Direction& x)
              { return f(x,s) + gDeriv(x,s) - sigma*g(x,s); },
              grids[i]->leafGridView());
        auto gApprox = Functions::makeDiscreteGlobalBasisFunction<double>(
              feBasisInterior, boundaryValues[i]);

        auto localGExact = localFunction(gExact);
        auto localGApprox = localFunction(gApprox);
        auto localView = feBasisInterior.localView();
        auto localIndexSet = feBasisInterior.localIndexSet();

        double rhsError_i = 0.;
        std::vector<std::tuple<EntitySeed, double>> errorEstimates;
        errorEstimates.reserve(feBasisInterior->gridView().size(0));
        for(const auto& e : elements(feBasisInterior->gridView())) {
          localView.bind(e);
          localIndexSet.bind(localView);
          localGExact.bind(e);
          localGApprox.bind(e);

          size_t quadratureOrder = 2*localView.tree().finiteElement()
                                              .localBasis().order()  + 4;
          const Dune::QuadratureRule<double, dim>& quad =
                Dune::QuadratureRules<double, dim>::rule(e.type(),
                                                         quadratureOrder);
          auto geometry = e.geometry();
          double local_error = 0.;
          for (size_t pt=0, qsize=quad.size(); pt < qsize; pt++) {
            const FieldVector<double,dim>& quadPos = quad[pt].position();
            const double diff = localGExact(quadPos) - localGApprox(quadPos);
            local_error += diff*diff
                         * geometry.integrationElement(quadPos)
                         * quad[pt].weight();
          }
          errorEstimates.emplace_back(e.seed(), local_error);
          rhsError_i += local_error;
        }
        rhsError_i = std::sqrt(rhsError_i);

        if(rhsError_i <= accuracy / numS ||
            refinement == maxNumberOfRefinements) {
          break;
        } else {
          ErrorTools::DoerflerMarking(*grids[i], 0.2,
              std::move(errorEstimates));
          auto oldGridData = attachDataToGrid(feBasisInterior,
                                              rhsFunctional[i]);
          grids[i]->preAdapt();
          grids[i]->adapt();
          grids[i]->postAdapt();
          detail::updateSpaces(*solutionSpaces[i], grids[i]->leafGridView());
          // TODO: also need to update test spaces.
          auto newGridData = restoreDataToRefinedGrid(feBasisInterior,
                                                      oldGridData);
          rhsFunctional[i].resize(newGridData.size(),
                                  false /* don't copy old values */);
          for(size_t j = 0, jmax = newGridData.size(); j < jmax; j++) {
            rhsFunctional[i][j] = newGridData[j];
          }
        }
      }

      // Add boundaryValues[i] to first feBasisInterior.size() entries of
      // rhsFunctional[i].
      using Iterator = std::decay_t<decltype(rhsFunctional[i].begin())>;
      for(Iterator rIt=rhsFunctional[i].begin(),
                   rEnd=rhsFunctional[i].begin()+feBasisInterior.size(),
                   gIt=boundaryValues[i].begin();
          rIt!=rEnd; ++rIt, ++gIt) {
        *rIt += *gIt;
      }
    }
  }
} // end namespace detail
#endif

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class HostGrid, class F, class G, class GDeriv, class Kernel>
void Periter<ScatteringKernelApproximation, RHSApproximation>::solve(
           HostGrid& hostGrid,
           const F& f,
           const G& g,
           const GDeriv& gDeriv,
           double sigma,
           const Kernel& kernel,
           unsigned int numS,
           double rho,
           double CT,
           double targetAccuracy,
           unsigned int maxNumberOfIterations,
           PlotSolutions plotSolutions) {
  if(plotSolutions == PlotSolutions::plotLastIteration) {
    std::cerr
        << "Plotting of only the last iteration is not implemented yet!\n";
    std::abort();
  }
  static_assert(std::is_same<RHSApproximation, FeRHSandBoundary>::value
      || std::is_same<RHSApproximation, ApproximateRHSandBoundary>::value,
      "Unknown type provided for RHSApproximation!\n"
      "Should be either FeRHSandBoundary or ApproximateRHSandBoundary.");
  constexpr bool rhsIsFeFunction
      = std::is_same<RHSApproximation, FeRHSandBoundary>::value;

  const unsigned int dim = HostGrid::dimension;
  using Grid = SubGrid<dim, HostGrid, false>;
  using LeafGridView = typename Grid::LeafGridView;
  using HostGridView = typename HostGrid::LeafGridView;
  using Geometry = typename LeafGridView::template Codim<0>::Geometry;
  using Domain = typename Geometry::GlobalCoordinate;
  using Direction = FieldVector<double, dim>;
  using GridIdSet = std::set<typename HostGrid::GlobalIdSet::IdType>;

  ///////////////////////////////////
  // To print information
  ///////////////////////////////////
  std::ofstream ofs("output_rad_trans");

  ////////////////////////////////////////////
  // Handle directions of discrete ordinates
  ////////////////////////////////////////////
  // Vector of directions: sVector
  std::vector< Direction > sVector(numS);
  for(unsigned int i = 0; i < numS; ++i)
  {
    using namespace boost::math::constants;
    sVector[i] = {cos(2*pi<double>()*i/numS),
                  sin(2*pi<double>()*i/numS)};
  }

  std::vector<std::unique_ptr<Grid>> grids;
  grids.reserve(sVector.size());
  for(size_t i = 0, imax = sVector.size(); i < imax; i++) {
    std::unique_ptr<Grid> gr = std::make_unique<Grid>(hostGrid);
    gr->createBegin();
    gr->insertLevel(hostGrid.maxLevel());
    gr->createEnd();
    gr->setMaxLevelDifference(1);
    grids.emplace_back(std::move(gr));
  }

  std::vector<GridIdSet> gridIdSets;
  gridIdSets.reserve(grids.size());
  for(size_t i = 0, imax = grids.size(); i < imax; i++) {
    gridIdSets.push_back(saveSubGridToIdSet(*grids[i]));
  }

  //////////////////////////////////////////////////////////////////////
  //   Choose finite element spaces for the solution and test functions
  //////////////////////////////////////////////////////////////////////

  typedef Functions::LagrangeDGBasis<LeafGridView, 1> FEBasisInterior; // u
  typedef Functions::HangingNodeP2NodalBasis<LeafGridView> FEBasisTrace; // u^

  using FEBasisTest = Functions::LagrangeDGBasis<LeafGridView, 4>;
  using FEBasisTestEnriched = FEBasisTest;

  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////
  typedef BlockVector<FieldVector<double,1> > VectorType;
  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

  /////////////////////////////////////////////////
  //   Solution vectors
  /////////////////////////////////////////////////
  std::vector<VectorType> x(numS);
  std::vector<decltype(
      make_space_tuple<FEBasisInterior, FEBasisTrace>
          (std::declval<LeafGridView>()))> solutionSpaces;
  solutionSpaces.reserve(numS);
  std::vector<decltype(
      make_space_tuple<FEBasisTest>(std::declval<LeafGridView>()))> testSpaces;
  testSpaces.reserve(numS);
  std::vector<decltype(
      make_space_tuple<FEBasisTestEnriched>(std::declval<LeafGridView>()))>
        testSpacesEnriched;
  testSpacesEnriched.reserve(numS);

  for(unsigned int i = 0; i < numS; ++i)
  {
    auto gridView = grids[i]->leafGridView();
    solutionSpaces.push_back(
        make_space_tuple<FEBasisInterior, FEBasisTrace>(gridView));
    // v enriched
    testSpaces.push_back(make_space_tuple<FEBasisTest>(gridView));
    testSpacesEnriched.push_back(
        make_space_tuple<FEBasisTestEnriched>(gridView));

    x[i].resize(std::get<0>(*solutionSpaces[i]).size()
               + std::get<1>(*solutionSpaces[i]).size());
    x[i] = 0;
  }

  ScatteringKernelApproximation kernelApproximation(kernel, numS);

  /////////////////////////////////////////////////////////
  //  Fixed-point iterations
  /////////////////////////////////////////////////////////
  // TODO: A priori estimate for the accuracy of our solution:
  double accuracy = 1.;
  // η_n:
  double eta = 1;
  // TODO: estimate norm of rhs f
  const double fnorm = 1;
  // ρ̄:
  const double rhobar = (1./rho > 2*rho)? (1./rho) : (2*rho);

  // CT*kappa1 + (1+CT)*kappa2 + 2*kappa3 = 1.
  const double kappa1 = rhsIsFeFunction? 1./(2.*CT) : 1./(3.*CT);
  const double kappa2 = rhsIsFeFunction? 0.         : 1./(3.*(1+CT));
  const double kappa3 = rhsIsFeFunction? 1./4.      : 1./6.;

  ofs << "Periter with " << numS
      << " directions, rho = " << rho << ", CT = " << CT
      << ", kappa1 = " << kappa1
      << ", kappa2 = " << kappa2
      << ", kappa3 = " << kappa3
      << std::endl;

  for(unsigned int n = 0; accuracy > targetAccuracy
                          && n < maxNumberOfIterations; ++n)
  {
    ofs << "\nIteration " << n << std::endl;
    std::cout << "\nIteration " << n << std::endl << std::endl;

    for(size_t i = 0; i < numS; ++i) {
      grids[i] = restoreSubGridFromIdSet<Grid>(gridIdSets[i],
                                               hostGrid);
      detail::updateSpaces(*solutionSpaces[i], grids[i]->leafGridView());
      detail::updateSpaces(*testSpaces[i], grids[i]->leafGridView());
      detail::updateSpaces(*testSpacesEnriched[i], grids[i]->leafGridView());
    }

    std::chrono::steady_clock::time_point startScatteringApproximation
        = std::chrono::steady_clock::now();

    double kappaNorm = 1.;
    double uNorm = 0.;
    for(size_t i=0; i<numS; ++i) {
      const double uiNorm =
        ErrorTools::l2norm(std::get<FEBasisInterior>(*solutionSpaces[i]), x[i]);
      uNorm += uiNorm * uiNorm;
    }
    uNorm = std::sqrt(uNorm / numS);

    using FEBasisHostInterior
        = changeGridView_t<FEBasisInterior, HostGridView>;
    using RHSData = std::decay_t<decltype(attachDataToSubGrid(
              std::declval<FEBasisTest>(),
              std::declval<FEBasisHostInterior>(),
              std::declval<VectorType>()))>;
    std::vector<RHSData> rhsData;
    rhsData.reserve(numS);
    {
      FEBasisHostInterior hostGridGlobalBasis(hostGrid.leafGridView());
      std::vector<VectorType> rhsFunctional =
          apply_scattering (
            kernelApproximation, x, solutionSpaces, hostGridGlobalBasis,
            kappa1 * eta / (kappaNorm * uNorm));

      // TODO: restore grids from gridIdSets and update spaces
      //       shouldn't be necessary, as long as RHSApproximation ==
      //       FERHSandBoundary
      detail::approximate_rhs_and_bv (
          rhsFunctional,
          grids,
          kappa2*eta,
          hostGridGlobalBasis,
          sVector,
          f, g, gDeriv, sigma, RHSApproximation{});
      for(size_t i = 0; i < numS; i++) {
        FEBasisTest& subGridGlobalBasis
            = std::get<FEBasisTest>(*testSpaces[i]);
        rhsData.push_back(
          attachDataToSubGrid(
            subGridGlobalBasis,
            hostGridGlobalBasis,
            rhsFunctional[i]));
      }
    }

    std::chrono::steady_clock::time_point endScatteringApproximation
        = std::chrono::steady_clock::now();


    ////////////////////////////////////////////////////
    // Inner loop
    ////////////////////////////////////////////////////
    double accumulatedAPosterioriError = 0.;
    for(unsigned int i = 0; i < numS; ++i)
    {
      const Direction s = sVector[i];

      grids[i] = restoreSubGridFromIdSet<Grid>(gridIdSets[i],
                                               hostGrid);
      detail::updateSpaces(*solutionSpaces[i], grids[i]->leafGridView());
      detail::updateSpaces(*testSpaces[i], grids[i]->leafGridView());
      detail::updateSpaces(*testSpacesEnriched[i], grids[i]->leafGridView());

      auto bilinearForm =
        make_BilinearForm(testSpaces[i], solutionSpaces[i],
            make_tuple(
                make_IntegralTerm<0,0,IntegrationType::valueValue,
                                      DomainOfIntegration::interior>(sigma),
                make_IntegralTerm<0,0,IntegrationType::gradValue,
                                      DomainOfIntegration::interior>(-1., s),
                make_IntegralTerm<0,1,IntegrationType::normalVector,
                                      DomainOfIntegration::face>(1., s)));
      auto bilinearFormEnriched =
          replaceTestSpaces(bilinearForm, testSpacesEnriched[i]);
      auto innerProduct =
        make_InnerProduct(testSpaces[i],
            make_tuple(
                make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                      DomainOfIntegration::interior>(1., s),
                make_IntegralTerm<0,0,IntegrationType::travelDistanceWeighted,
                                      DomainOfIntegration::face>(1., s)));
      auto innerProductEnriched =
          replaceTestSpaces(innerProduct, testSpacesEnriched[i]);


      typedef GeometryBuffer<typename LeafGridView::template Codim<0>::Geometry>
          GeometryBuffer_t;

      GeometryBuffer_t geometryBuffer;
      auto systemAssembler =
          make_DPGSystemAssembler(bilinearForm,
                                  innerProduct,
                                  geometryBuffer);

      const unsigned int maxNumberOfInnerIterations = 3;
      double aposterioriErr;
      for(unsigned int nRefinement = 0; ; )
        // At the end of the loop, we will break if
        // aposterioriErr < kapp3*eta TODO: / numS ?
        // or ++nRefinement >= maxNumberOfInnerIterations
        // thus the inner loop terminates eventually.
      {
        std::cout << "Direction " << i
                  << ", inner iteration " << nRefinement << '\n';

        // Determine Dirichlet dofs for u^ (inflow boundary)
        std::vector<bool> dirichletNodesInflow;
        // Contribution of inflow boundary for the rhs
        std::vector<double> rhsInflowContrib;
        {
          BoundaryTools::getInflowBoundaryMask(std::get<1>(*solutionSpaces[i]),
                                                dirichletNodesInflow,
                                                s);

          auto gSfixed = [s](const Domain& x){ return 0.;};
          BoundaryTools::getBoundaryValue(std::get<1>(*solutionSpaces[i]),
                                          rhsInflowContrib,
                                          gSfixed);
        }

        VectorType rhsFunctional;
        {
          FEBasisTest& feBasisTest = std::get<FEBasisTest>(*testSpaces[i]);
          auto newGridData
              = rhsData[i].restoreDataToRefinedSubGrid(feBasisTest);
          rhsFunctional.resize(newGridData.size(),
                               false /* don't copy old values */);
          for(size_t j = 0, jmax = newGridData.size(); j < jmax; j++) {
            rhsFunctional[j] = newGridData[j];
          }
        }

        VectorType rhs;
        MatrixType stiffnessMatrix;

        /////////////////////////////////////////////////////////
        //  Assemble the systems
        /////////////////////////////////////////////////////////
        // loop of the discrete ordinates
        {
          auto rhsFunction = make_LinearForm(
                systemAssembler.getTestSearchSpaces(),
                std::make_tuple(
                  make_LinearFunctionalTerm<0, DomainOfIntegration::interior>
                    (rhsFunctional, std::get<0>(*testSpaces[i]))));
          systemAssembler.assembleSystem(
              stiffnessMatrix, rhs,
              rhsFunction);
          systemAssembler.template applyDirichletBoundary<1>
              (stiffnessMatrix,
              rhs,
              dirichletNodesInflow,
              rhsInflowContrib);
#if 1
          systemAssembler.template defineCharacteristicFaces<1,dim>(
              stiffnessMatrix,
              rhs, s);
#endif
        }

        ////////////////////////////////////
        //   Initialize solution vector
        ////////////////////////////////////
        x[i].resize(std::get<0>(*solutionSpaces[i]).size()
                    + std::get<1>(*solutionSpaces[i]).size());
        x[i] = 0;

        ////////////////////////////
        //   Compute solution
        ////////////////////////////

        std::cout << "rhs size = " << rhs.size()
                  << " matrix size = " << stiffnessMatrix.N()
                              << " x " << stiffnessMatrix.M()
                  << " solution size = " << x[i].size() << std::endl;


        const int verbosity = 0; // 0: not verbose; >0: verbose
        UMFPack<MatrixType> umfPack(stiffnessMatrix, verbosity);
        InverseOperatorResult statistics;
        umfPack.apply(x[i], rhs, statistics);

        ////////////////////////////////////
        //  A posteriori error
        ////////////////////////////////////
        std::cout << "Compute a posteriori error\n";

        // We compute the a posteriori error
        // - We compute the rhs with the enriched test space ("rhs=f(v)")
        // -- Contribution of the source term f that has an analytic expression
        // -- Contribution of the scattering term
        {
          auto rhsAssemblerEnriched = make_RhsAssembler(testSpacesEnriched[i]);
          auto rhsFunction = make_LinearForm(
                rhsAssemblerEnriched.getTestSpaces(),
                std::make_tuple(
                  make_LinearFunctionalTerm<0, DomainOfIntegration::interior>
                    (rhsFunctional,
                    std::get<0>(*systemAssembler.getTestSearchSpaces()))));
          rhsAssemblerEnriched.assembleRhs(rhs, rhsFunction);
        }
        // - Computation of the a posteriori error
        using EntitySeed = typename Grid::template Codim<0>::Entity::EntitySeed;
        std::vector<std::tuple<EntitySeed, double>> errorEstimates
            = ErrorTools::residual(bilinearFormEnriched,
                                   innerProductEnriched,
                                   x[i], rhs);
        aposterioriErr = std::sqrt(std::accumulate(
              errorEstimates.cbegin(), errorEstimates.cend(), 0.,
              [](double acc, const std::tuple<EntitySeed, double>& t)
              {
                return acc + std::get<1>(t);
              }));

        static_assert(!is_RefinedFiniteElement<FEBasisInterior>::value,
            "Functions::interpolate won't work for refined finite elements");
        {
          auto& feBasisInterior = std::get<0>(*solutionSpaces[i]);
          VectorType gInterpolation(feBasisInterior.size());
          Functions::interpolate(feBasisInterior, gInterpolation,
              [&g,s](const Direction& x) { return g(x,s); });

          // Add gInterpolation to first feBasisInterior.size() entries of x.
          using Iterator = std::decay_t<decltype(x[i].begin())>;
          for(Iterator xIt=x[i].begin(),
                       xEnd=x[i].begin()+feBasisInterior.size(),
                       gIt=gInterpolation.begin(); xIt!=xEnd; ++xIt, ++gIt) {
            *xIt += *gIt;
          }
          // TODO: Add (interpolation of) g to theta part of x[i]?
        }

        ofs << "Iteration " << n << '.' << nRefinement
            << "for direction " << i << ": "
            << "A posteriori estimation of || (u,trace u) - (u_fem,theta) || = "
            << aposterioriErr << ", grid level: " << grids[i]->maxLevel()
            << ", number of DOFs: " << x[i].size()
            << ", applying the kernel took "
            << std::chrono::duration_cast<std::chrono::microseconds>
                (endScatteringApproximation - startScatteringApproximation)
                .count()
            << "us, " << kernelApproximation.info()
            << '\n';
        std::cout << '\n';

        std::cout << "\nStatistics at end of inner iteration:\n";
        std::cout << "Grid level: " << grids[i]->maxLevel() << '\n';
        std::cout << "A posteriori error: " << aposterioriErr << '\n';

        if(++nRefinement >= maxNumberOfInnerIterations
            || aposterioriErr <= kappa3*eta) {
          gridIdSets[i] = saveSubGridToIdSet(*grids[i]);
          break;
        } else {
          const double ratio = .2;
          ErrorTools::DoerflerMarking(*grids[i], ratio,
                                      std::move(errorEstimates));
          grids[i]->preAdapt();
          grids[i]->adapt();
          grids[i]->postAdapt();
          detail::updateSpaces(*solutionSpaces[i], grids[i]->leafGridView());
          detail::updateSpaces(*testSpaces[i], grids[i]->leafGridView());
          detail::updateSpaces(*testSpacesEnriched[i], grids[i]->leafGridView());
        }
      }
      accumulatedAPosterioriError += aposterioriErr * aposterioriErr;
    }
    accumulatedAPosterioriError
        = std::sqrt(accumulatedAPosterioriError / static_cast<double>(numS));
    const size_t accumulatedDoFs = std::accumulate(x.cbegin(), x.cend(),
        static_cast<size_t>(0),
        [](size_t acc, auto vec) { return acc + vec.size(); });
    ofs << "Error at end of Iteration " << n << ": "
        << accumulatedAPosterioriError << ", using "
        << accumulatedDoFs << " DoFs, applying the kernel took "
        << std::chrono::duration_cast<std::chrono::microseconds>
            (endScatteringApproximation - startScatteringApproximation)
            .count()
        << "us, " << kernelApproximation.info()
        << '\n';
    std::cout << "Error at end of Iteration " << n << ": "
              << accumulatedAPosterioriError << ", using "
              << accumulatedDoFs << " DoFs\n";

    accuracy = std::pow(rho, n) * CT * fnorm + 2*eta;
    eta /= rhobar;

    if(plotSolutions == PlotSolutions::plotOuterIterations) {
      ////////////////////////////////////////////////////////////////////////
      //  Write result to VTK file
      //  We need to subsample, because VTK cannot natively display
      //  real second-order functions
      ////////////////////////////////////////////////////////////////////////
      std::cout << "Print solutions:\n";

      for(unsigned int i = 0; i < numS; ++i)
      {
        grids[i] = restoreSubGridFromIdSet<Grid>(gridIdSets[i],
                                                 hostGrid);
        detail::updateSpaces(*solutionSpaces[i], grids[i]->leafGridView());
        detail::updateSpaces(*testSpaces[i], grids[i]->leafGridView());
        detail::updateSpaces(*testSpacesEnriched[i], grids[i]->leafGridView());

        FEBasisInterior& feBasisInterior
            = std::get<FEBasisInterior>(*solutionSpaces[i]);
        FEBasisTrace& feBasisTrace
            = std::get<FEBasisTrace>(*solutionSpaces[i]);

        std::cout << "Direction " << i << '\n';

        std::string name = std::string("u_rad_trans_n")
                        + std::to_string(n)
                        + std::string("_s")
                        + std::to_string(i);
        FunctionPlotter uPlotter(name);
        uPlotter.plot("u", x[i], feBasisInterior, 0, 0);
        name = std::string("theta_rad_trans_n")
                        + std::to_string(n)
                        + std::string("_s")
                        + std::to_string(i);
        FunctionPlotter thetaPlotter(name);
        thetaPlotter.plot("theta", x[i], feBasisTrace, 2,
                          feBasisInterior.size());
      }
    }
  }
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class SolutionSpaces, class HostGridBasis>
std::vector<typename Periter<ScatteringKernelApproximation,
                             RHSApproximation>::VectorType>
Periter<ScatteringKernelApproximation, RHSApproximation>::apply_scattering(
      ScatteringKernelApproximation& kernelApproximation,
      const std::vector<VectorType>& x,
      const std::vector<std::shared_ptr<SolutionSpaces>>& solutionSpaces,
      const HostGridBasis& hostGridBasis,
      double accuracy) {
  kernelApproximation.setAccuracy(accuracy);

  using FEBasisInterior = std::tuple_element_t<0, SolutionSpaces>;

  const size_t numS = x.size();
  // Interpolate x[i] to hostGridBasis.
  std::vector<VectorType> xHost(numS);
  for(size_t i = 0; i < numS; ++i) {
    FEBasisInterior& feBasisInterior = std::get<0>(*solutionSpaces[i]);
    interpolateFromSubGrid(
        feBasisInterior, x[i],
        hostGridBasis, xHost[i]);
  }

  auto scatteringAssembler =
      make_ApproximateScatteringAssembler(hostGridBasis,
                                          kernelApproximation);
  std::vector<VectorType> rhsFunctional(numS);
  for(size_t i = 0; i < numS; ++i) {
    scatteringAssembler.precomputeScattering(
        rhsFunctional[i],
        xHost, i);
  }

  return rhsFunctional;
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH
