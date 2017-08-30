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
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

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
    const std::set<typename SubGrid::HostGridType::GlobalIdSet::IdType>& idSet,
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

  using Direction = typename ScatteringKernelApproximation::Direction;

  /**
   * Solve a radiative transfer problem using the Periter algorithm
   *
   * \param hostGrid
   * \param f  right hand side function
   * \param g  lifting of the boundary values
   * \param gDeriv  derivative of g in direction s
   * \param sigma   absorbtion coefficient
   * \param kernel  the scattering kernel, e.g. a Henyey–Greenstein kernel
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
  template<class SolutionSpaces, class HostGrid, class Grid>
  static std::vector<VectorType> apply_scattering(
      ScatteringKernelApproximation& kernelApproximation,
      const std::vector<VectorType>& x,
      const std::vector<std::shared_ptr<SolutionSpaces>>& solutionSpaces,
      const HostGrid& hostGrid,
      std::vector<Direction>& sVector,
      std::vector<std::unique_ptr<Grid>>& grids,
      double accuracy);

  /**
   * create grids for new directions created in apply_scattering
   */
  template<class Grid>
  static void create_new_grids(
      std::vector<std::unique_ptr<Grid>>& grids,
      size_t numNewGrids);

  /**
   * insert new gridIdSets after apply_scattering added new grids
   */
  template<class GridIdSet, class Grid>
  static void create_new_gridIdSets(
      std::vector<GridIdSet>& gridIdSets,
      const std::vector<std::unique_ptr<Grid>>& grids);

  /**
   * insert new spaces after apply_scattering added new grids
   */
  template<class Grid, class... Spaces>
  static void create_new_spaces(
      std::vector<std::shared_ptr<std::tuple<Spaces...>>>& spaces,
      const std::vector<std::unique_ptr<Grid>>& grids);
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
    const size_t numS = sVector.size();
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

  ///////////////////////////////////
  // Parameters for adaptivity
  ///////////////////////////////////

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

  ////////////////////////////////////////////
  // Handle directions of discrete ordinates
  ////////////////////////////////////////////

  // TODO: The accuracy also depends on the kappa from K = κG and on \|u\|.
  //       Adding a factor 1/4. to compensate for that.
  ScatteringKernelApproximation kernelApproximation(kernel,
                                                    kappa1*targetAccuracy/4.);

  ofs << "Periter with up to " << kernelApproximation.maxNumS()
      << " directions, rho = " << rho << ", CT = " << CT
      << ", kappa1 = " << kappa1
      << ", kappa2 = " << kappa2
      << ", kappa3 = " << kappa3
      << '\n';

  // As the solution u we use for the initial scattering is 0, and the
  // formula for the accuracy contains a 1/\|u\|, we set the initial
  // accuracy to a large enough value.
  std::vector<Direction> sVector(kernelApproximation.setAccuracy(1e5));
  unsigned int numS = sVector.size();

  std::vector<std::unique_ptr<Grid>> grids;
  grids.reserve(sVector.size()/kernelApproximation.numSperInterval);
  for(size_t i = 0, imax = sVector.size()/kernelApproximation.numSperInterval;
      i < imax; i++) {
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

  using FEBasisTest = Functions::PQkDGRefinedDGBasis<LeafGridView, 1, 3>;
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
  solutionSpaces.reserve(grids.size());
  std::vector<decltype(
      make_space_tuple<FEBasisTest>(std::declval<LeafGridView>()))> testSpaces;
  testSpaces.reserve(grids.size());
  std::vector<decltype(
      make_space_tuple<FEBasisTestEnriched>(std::declval<LeafGridView>()))>
        testSpacesEnriched;
  testSpacesEnriched.reserve(grids.size());

  for(unsigned int i = 0; i < grids.size(); ++i)
  {
    auto gridView = grids[i]->leafGridView();
    solutionSpaces.push_back(
        make_space_tuple<FEBasisInterior, FEBasisTrace>(gridView));
    // v enriched
    testSpaces.push_back(make_space_tuple<FEBasisTest>(gridView));
    testSpacesEnriched.push_back(
        make_space_tuple<FEBasisTestEnriched>(gridView));

    for(unsigned int j = i*kernelApproximation.numSperInterval,
                     jmax = (i+1)*kernelApproximation.numSperInterval;
        j < jmax; j++) {
      x[j].resize(std::get<0>(*solutionSpaces[i]).size()
                  + std::get<1>(*solutionSpaces[i]).size());
      x[j] = 0;
    }
  }

  /////////////////////////////////////////////////////////
  //  Fixed-point iterations
  /////////////////////////////////////////////////////////

  // TODO: A priori estimate for the accuracy of our solution:
  double accuracy = 1.;

  for(unsigned int n = 0; accuracy > targetAccuracy
                          && n < maxNumberOfIterations; ++n)
  {
    ofs << "\nIteration " << n << std::endl;
    std::cout << "\nIteration " << n << std::endl << std::endl;

    for(size_t i = 0; i < grids.size(); ++i) {
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
    for(size_t i=0; i<solutionSpaces.size(); ++i) {
      for(unsigned int j = i*kernelApproximation.numSperInterval,
                      jmax = (i+1)*kernelApproximation.numSperInterval;
          j < jmax; j++) {
        const double ujNorm =
          ErrorTools::l2norm(std::get<FEBasisInterior>(*solutionSpaces[i]),
                             x[j]);
        uNorm += ujNorm * ujNorm;
      }
    }
    uNorm = std::sqrt(uNorm / numS);
    // To prevent division by zero.
    if(uNorm == 0.) uNorm = 1e-5;

    using FEBasisHostInterior
        = changeGridView_t<FEBasisInterior, HostGridView>;
    using RHSData = std::decay_t<decltype(attachDataToSubGrid(
              std::declval<FEBasisTest>(),
              std::declval<FEBasisHostInterior>(),
              std::declval<VectorType>()))>;
    std::vector<RHSData> rhsData;
    {
      FEBasisHostInterior hostGridGlobalBasis(hostGrid.leafGridView());
      std::vector<VectorType> rhsFunctional =
          apply_scattering (
            kernelApproximation, x, solutionSpaces, hostGridGlobalBasis,
            sVector, grids,
            kappa1 * eta / (kappaNorm * uNorm));
      create_new_gridIdSets(gridIdSets, grids);
      create_new_spaces(solutionSpaces, grids);
      create_new_spaces(testSpaces, grids);
      create_new_spaces(testSpacesEnriched, grids);

      numS = sVector.size();
      x.resize(numS);
      rhsData.reserve(numS);

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
      size_t i = 0;
      for(auto& testSpace : testSpaces) {
        const FEBasisTest& subGridGlobalBasis
            = std::get<FEBasisTest>(*testSpace);
        for(size_t imax = i + kernelApproximation.numSperInterval;
            i < imax; i++) {
          rhsData.push_back(
            attachDataToSubGrid(
              subGridGlobalBasis,
              hostGridGlobalBasis,
              rhsFunctional[i]));
        }
      }
    }

    std::chrono::steady_clock::time_point endScatteringApproximation
        = std::chrono::steady_clock::now();


    ////////////////////////////////////////////////////
    // Inner loop
    ////////////////////////////////////////////////////
    double accumulatedAPosterioriError = 0.;
    for(unsigned int i = 0; i < grids.size(); ++i)
    {
      grids[i] = restoreSubGridFromIdSet<Grid>(gridIdSets[i],
                                               hostGrid);
      detail::updateSpaces(*solutionSpaces[i], grids[i]->leafGridView());
      detail::updateSpaces(*testSpaces[i], grids[i]->leafGridView());
      detail::updateSpaces(*testSpacesEnriched[i], grids[i]->leafGridView());

      const unsigned int maxNumberOfInnerIterations = 16;
      double aposterioriErr;
      unsigned int nRefinement = 0;
      // TODO: refine grid for all all directions in same interval at once.
      for( ; ; )
        // At the end of the loop, we will break if
        // aposterioriErr < kapp3*eta (pointwise in s)
        // or ++nRefinement >= maxNumberOfInnerIterations
        // thus the inner loop terminates eventually.
      {
        aposterioriErr = 0.;

        for(unsigned int j = i*kernelApproximation.numSperInterval,
                        jmax = (i+1)*kernelApproximation.numSperInterval;
            j < jmax; j++)
        {
          std::cout << "Direction " << j
                    << ", inner iteration " << nRefinement << '\n';

          const Direction s = sVector[j];

          auto bilinearForm =
            make_BilinearForm(testSpaces[i], solutionSpaces[i],
                make_tuple(
                    make_IntegralTerm<0,0,IntegrationType::valueValue,
                                          DomainOfIntegration::interior>(sigma),
                    make_IntegralTerm<0,0,
                                      IntegrationType::gradValue,
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
                    make_IntegralTerm<0,0,
                                      IntegrationType::travelDistanceWeighted,
                                      DomainOfIntegration::face>(1., s)));
          auto innerProductEnriched =
              replaceTestSpaces(innerProduct, testSpacesEnriched[i]);


          using Geometry = typename LeafGridView::template Codim<0>::Geometry;
          using GeometryBuffer_t = GeometryBuffer<Geometry>;

          GeometryBuffer_t geometryBuffer;
          auto systemAssembler =
              make_DPGSystemAssembler(bilinearForm,
                                      innerProduct,
                                      geometryBuffer);

          // Determine Dirichlet dofs for u^ (inflow boundary)
          std::vector<bool> dirichletNodesInflow;
          // Contribution of inflow boundary for the rhs
          std::vector<double> rhsInflowContrib;
          {
            BoundaryTools::getInflowBoundaryMask(
                        std::get<1>(*solutionSpaces[i]),
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
                = rhsData[j].restoreDataToRefinedSubGrid(feBasisTest);
            rhsFunctional.resize(newGridData.size(),
                                 false /* don't copy old values */);
            for(size_t k = 0, kmax = newGridData.size(); k < kmax; k++) {
              rhsFunctional[k] = newGridData[k];
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
#if 0
            systemAssembler.template defineCharacteristicFaces<1>(
                stiffnessMatrix,
                rhs, s);
#endif
          }

          ////////////////////////////////////
          //   Initialize solution vector
          ////////////////////////////////////
          x[j].resize(std::get<0>(*solutionSpaces[i]).size()
                      + std::get<1>(*solutionSpaces[i]).size());
          x[j] = 0;

          ////////////////////////////
          //   Compute solution
          ////////////////////////////

          std::cout << "rhs size = " << rhs.size()
                    << " matrix size = " << stiffnessMatrix.N()
                                << " x " << stiffnessMatrix.M()
                    << " solution size = " << x[j].size() << std::endl;


          const int verbosity = 0; // 0: not verbose; >0: verbose
          UMFPack<MatrixType> umfPack(stiffnessMatrix, verbosity);
          InverseOperatorResult statistics;
          umfPack.apply(x[j], rhs, statistics);

          ////////////////////////////////////
          //  A posteriori error
          ////////////////////////////////////
          std::cout << "Compute a posteriori error\n";

          // We compute the a posteriori error
          // - We compute the rhs with the enriched test space ("rhs=f(v)")
          // -- Contribution of the source term f that has an
          //    analytic expression
          // -- Contribution of the scattering term
          {
            auto rhsAssemblerEnriched
              = make_RhsAssembler(testSpacesEnriched[i]);
            auto rhsFunction = make_LinearForm(
                  rhsAssemblerEnriched.getTestSpaces(),
                  std::make_tuple(
                    make_LinearFunctionalTerm<0, DomainOfIntegration::interior>
                      (rhsFunctional,
                      std::get<0>(*systemAssembler.getTestSearchSpaces()))));
            rhsAssemblerEnriched.assembleRhs(rhs, rhsFunction);
          }
          // - Computation of the a posteriori error
          using EntitySeed
              = typename Grid::template Codim<0>::Entity::EntitySeed;
          std::vector<std::tuple<EntitySeed, double>> errorEstimates
              = ErrorTools::residual(bilinearFormEnriched,
                                    innerProductEnriched,
                                    x[j], rhs);
          aposterioriErr = std::sqrt(std::accumulate(
                errorEstimates.cbegin(), errorEstimates.cend(),
                aposterioriErr * aposterioriErr,
                [](double acc, const std::tuple<EntitySeed, double>& t)
                {
                  return acc + std::get<1>(t);
                }));
          const double ratio = .6;
          ErrorTools::DoerflerMarking(*grids[i], ratio,
                                      std::move(errorEstimates));

          static_assert(!is_RefinedFiniteElement<FEBasisInterior>::value,
              "Functions::interpolate won't work for refined finite elements");
          {
            auto& feBasisInterior = std::get<0>(*solutionSpaces[i]);
            VectorType gInterpolation(feBasisInterior.size());
            Functions::interpolate(feBasisInterior, gInterpolation,
                [&g,s](const Direction& x) { return g(x,s); });

            // Add gInterpolation to first feBasisInterior.size() entries of x.
            using Iterator = std::decay_t<decltype(x[j].begin())>;
            for(Iterator xIt=x[j].begin(),
                        xEnd=x[j].begin()+feBasisInterior.size(),
                        gIt=gInterpolation.begin(); xIt!=xEnd; ++xIt, ++gIt) {
              *xIt += *gIt;
            }
            // TODO: Add (interpolation of) g to theta part of x[j]?
          }

          ofs << "Iteration " << n << '.' << nRefinement
              << " for direction " << j << ": "
              << "A posteriori estimation of || (u,trace u) - (u_fem,theta) || = "
              << aposterioriErr << ", grid level: " << grids[i]->maxLevel()
              << ", number of DOFs: " << x[j].size()
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
        }

        if(++nRefinement >= maxNumberOfInnerIterations
            || aposterioriErr <= kappa3*eta
                                 * kernelApproximation.numSperInterval) {
          gridIdSets[i] = saveSubGridToIdSet(*grids[i]);
          break;
        } else {
          grids[i]->preAdapt();
          grids[i]->adapt();
          grids[i]->postAdapt();
          detail::updateSpaces(*solutionSpaces[i], grids[i]->leafGridView());
          detail::updateSpaces(*testSpaces[i], grids[i]->leafGridView());
          detail::updateSpaces(*testSpacesEnriched[i],
                               grids[i]->leafGridView());
        }
      }
      accumulatedAPosterioriError += aposterioriErr * aposterioriErr;

      ofs << "after " << nRefinement << " transport solves, accuracy was "
          << ((aposterioriErr <=
               kappa3 * eta * kernelApproximation.numSperInterval)
              ?"reached.":"not reached.")
          << "\na posteriori error:   " << aposterioriErr
          << "\nprescribed tolerance: "
          << kappa3 * eta * kernelApproximation.numSperInterval
          << '\n';

      if(plotSolutions == PlotSolutions::plotOuterIterations) {
        //////////////////////////////////////////////////////////////////////
        //  Write result to VTK file
        //  We need to subsample, because VTK cannot natively display
        //  real second-order functions
        //////////////////////////////////////////////////////////////////////

        const FEBasisInterior& feBasisInterior
            = std::get<FEBasisInterior>(*solutionSpaces[i]);
        const FEBasisTrace& feBasisTrace
            = std::get<FEBasisTrace>(*solutionSpaces[i]);

        for(unsigned int j = i*kernelApproximation.numSperInterval,
                         jmax = (i+1)*kernelApproximation.numSperInterval;
            j < jmax; j++)
        {
          std::cout << "Plot solution for direction " << j << '\n';

          std::string name = std::string("u_rad_trans_n")
                            + std::to_string(n)
                            + std::string("_s")
                            // + std::to_string(j*stride);
                            + std::to_string(j);
          FunctionPlotter uPlotter(name);
          uPlotter.plot("u", x[j], feBasisInterior, 0, 0);
          name = std::string("theta_rad_trans_n")
                            + std::to_string(n)
                            + std::string("_s")
                            // + std::to_string(j*stride);
                            + std::to_string(j);
          FunctionPlotter thetaPlotter(name);
          thetaPlotter.plot("theta", x[j], feBasisTrace, 2,
                            feBasisInterior.size());
        }
      }
    }

    accumulatedAPosterioriError
        = std::sqrt(accumulatedAPosterioriError / static_cast<double>(numS));
    const size_t accumulatedDoFs = std::accumulate(x.cbegin(), x.cend(),
        static_cast<size_t>(0),
        [](size_t acc, auto vec) { return acc + vec.size(); });

    accuracy = std::pow(rho, n) * CT * fnorm + 2*eta;

    ofs << "Error at end of Iteration " << n << ": "
        << accumulatedAPosterioriError << ", using "
        << accumulatedDoFs << " DoFs, accuracy was " << accuracy
        << ", eta was " << eta
        << ", applying the kernel took "
        << std::chrono::duration_cast<std::chrono::microseconds>
            (endScatteringApproximation - startScatteringApproximation)
            .count()
        << "us, " << kernelApproximation.info()
        << '\n';
    std::cout << "Error at end of Iteration " << n << ": "
              << accumulatedAPosterioriError << ", using "
              << accumulatedDoFs << " DoFs\n";

    eta /= rhobar;
  }
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class SolutionSpaces, class HostGridBasis, class Grid>
std::vector<typename Periter<ScatteringKernelApproximation,
                             RHSApproximation>::VectorType>
Periter<ScatteringKernelApproximation, RHSApproximation>::apply_scattering(
      ScatteringKernelApproximation& kernelApproximation,
      const std::vector<VectorType>& x,
      const std::vector<std::shared_ptr<SolutionSpaces>>& solutionSpaces,
      const HostGridBasis& hostGridBasis,
      std::vector<Direction>& sVector,
      std::vector<std::unique_ptr<Grid>>& grids,
      double accuracy) {
  sVector = kernelApproximation.setAccuracy(accuracy);

  using FEBasisInterior = std::tuple_element_t<0, SolutionSpaces>;

  // Interpolate x[i] to hostGridBasis.
  std::vector<VectorType> xHost(x.size());
  {
    size_t i = 0;
    for(const auto& solutionSpace : solutionSpaces) {
      const FEBasisInterior& feBasisInterior = std::get<0>(*solutionSpace);
      for(size_t imax = i+kernelApproximation.numSperInterval;
          i < imax; ++i) {
        interpolateFromSubGrid(
            feBasisInterior, x[i],
            hostGridBasis, xHost[i]);
      }
    }
  }

  const auto scatteringAssembler =
      make_ApproximateScatteringAssembler(hostGridBasis,
                                          kernelApproximation);
  const size_t numS = sVector.size();
  std::vector<VectorType> rhsFunctional(numS);
  scatteringAssembler.precomputeScattering(rhsFunctional, xHost);

  create_new_grids(grids, numS/kernelApproximation.numSperInterval);

  return rhsFunctional;
}

/**
 * compute intersection of two SubGrids
 */
template<class SubGrid>
std::unique_ptr<SubGrid> intersectSubGrids(const SubGrid& subGrid1,
                                           const SubGrid& subGrid2)
{
  // assert(subGrid1.getHostGrid() == subGrid2.getHostGrid());
  std::unique_ptr<SubGrid> gr
      = std::make_unique<SubGrid>(subGrid1.getHostGrid());
  gr->createBegin();
  for(unsigned int level=0; level <= subGrid1.maxLevel(); ++level)
  {
    for (auto&& e : elements(subGrid1.levelGridView(level)))
    {
      const auto eHost = subGrid1.template getHostEntity<0>(e);
      if(subGrid2.template contains<0>(eHost))
      {
        gr->insertRaw(eHost);
      }
    }
  }
  gr->createEnd();
  gr->setMaxLevelDifference(1);
  return gr;
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class Grid>
void
Periter<ScatteringKernelApproximation, RHSApproximation>::
create_new_grids(
      std::vector<std::unique_ptr<Grid>>& grids,
      size_t numNewGrids)
{
  const size_t numOldGrids = grids.size();
  if(numOldGrids == numNewGrids) {
    // no new grids need to be added
    return;
  }
  std::vector<std::unique_ptr<Grid>> newGrids;
  newGrids.reserve(numNewGrids);
  const size_t numCopies = numNewGrids / numOldGrids;
  const size_t lastOldGrid = numOldGrids-1;
  for(size_t i = 0; i < lastOldGrid; i++) {
    newGrids.push_back(intersectSubGrids(*grids[i], *grids[i+1]));
    const Grid& intersectionGrid = *newGrids.back();
    for(size_t copies = 1; copies < numCopies; ++copies) {
      newGrids.push_back(copySubGrid(intersectionGrid));
    }
  }
  // last one needs to be handled extra, as one of the neighboring grids is
  // the first one (since the direction parameter lives on the unit sphere).
  newGrids.push_back(intersectSubGrids(*grids[lastOldGrid], *grids[0]));
  const Grid& intersectionGrid = *newGrids.back();
  for(size_t copies = 1; copies < numCopies; ++copies) {
    newGrids.push_back(copySubGrid(intersectionGrid));
  }

  std::swap(grids, newGrids);
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class GridIdSet, class Grid>
void
Periter<ScatteringKernelApproximation, RHSApproximation>::
create_new_gridIdSets(
      std::vector<GridIdSet>& gridIdSets,
      const std::vector<std::unique_ptr<Grid>>& grids)
{
  if(gridIdSets.size() != grids.size()) {
    gridIdSets.clear();
    gridIdSets.reserve(grids.size());
    for(auto grid = grids.cbegin(), end = grids.cend(); grid != end; ++grid) {
      gridIdSets.push_back(saveSubGridToIdSet(**grid));
    }
  }
}

template<class ScatteringKernelApproximation, class RHSApproximation>
template<class Grid, class... Spaces>
void
Periter<ScatteringKernelApproximation, RHSApproximation>::create_new_spaces(
      std::vector<std::shared_ptr<std::tuple<Spaces...>>>& spaces,
      const std::vector<std::unique_ptr<Grid>>& grids)
{
  if(spaces.size() != grids.size()) {
    spaces.clear();
    spaces.reserve(grids.size());
    for(auto grid = grids.cbegin(), end = grids.cend(); grid != end; ++grid) {
      spaces.push_back(make_space_tuple<Spaces...>((*grid)->leafGridView()));
    }
  }
}

} // end namespace Dune

#endif // DUNE_DPG_RADIATIVE_TRANSFER_PERITER_HH
