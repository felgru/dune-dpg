#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <cmath>
#include <iostream>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/dpg/functions/normedspaces.hh>
#include <dune/dpg/innerproduct.hh>
#include <dune/dpg/testspace_coefficient_matrix.hh>
#include <dune/functions/functionspacebases/bernsteindgrefineddgnodalbasis.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/istl/matrix.hh>

#include <Eigen/Core>
#include <Eigen/SVD>

#include <unistd.h>

using namespace Dune;

template<class InnerProduct, class TestLocalViews, class Element>
Eigen::MatrixXd stiffnessMatrix(InnerProduct& innerProduct,
                                TestLocalViews& testLocalViews,
                                const Element& e) {
  using namespace Dune::detail;
  bindLocalViews(testLocalViews, e);
  innerProduct.bind(testLocalViews);

  using MatrixType = Matrix<FieldMatrix<double,1,1>>;
  MatrixType sm;
  innerProduct.getLocalMatrix(sm);

  Eigen::MatrixXd stiffnessMatrix(sm.N(), sm.M());
  for(size_t i = 0; i < sm.N(); i++) {
    for(size_t j = 0; j < sm.N(); j++) {
      stiffnessMatrix(i, j) = sm[i][j][0];
    }
  }
  return stiffnessMatrix;
}

template<class InnerProduct>
double conditionOnFirstElement(InnerProduct& innerProduct)
{
  auto testSpaces = innerProduct.getTestSpaces();
  auto testLocalViews = Dune::detail::getLocalViews(*testSpaces);
  const auto gridView = std::get<0>(*testSpaces).gridView();
  for(const auto& e : elements(gridView)) {
    Eigen::MatrixXd ip
        = stiffnessMatrix(innerProduct, testLocalViews, e);
    Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::NoQRPreconditioner>
        svd(ip);

    // only visit first element
    return svd.singularValues()(0)
         / svd.singularValues()(ip.rows()-1);
  }
}

void printHelp(const char* name) {
  std::cerr << "Usage: " << name << " [-n] [-l m]\n"
               "Compute condition number of trial to test map for elements "
               "of decreasing size\n\n"
               " -n: normalize test functions\n"
               " -l m: maximal grid refinement level (default: 30)\n";
  exit(0);
}

int main(int argc, char *argv[]) {
  bool normalized = false;
  int numLevels = 30;
  int opt;
  while ((opt = getopt(argc,argv,"nl:h")) != EOF)
    switch(opt)
    {
      case 'n': normalized = true; break;
      case 'l': numLevels = atoi(optarg); break;
      default:
      case 'h':
      case '?':
        printHelp(argv[0]);
    }
  if(optind != argc) {
    printHelp(argv[0]);
  }

  constexpr int dim = 2;
  using Grid = UGGrid<dim>;

  const double theta = 1.;
  const FieldVector<double, dim> s = {std::cos(theta), std::sin(theta)};

  std::vector<std::shared_ptr<Grid>> grids;
  grids.reserve(numLevels);
  for(int level = 0; level < numLevels; level++) {
    const double size = std::exp2(-level);
    const FieldVector<double, dim> lower = {0, 0};
    const FieldVector<double, dim> upper = {size, size};
    const std::array<unsigned int, dim> numElements = {1, 1};

    const std::shared_ptr<Grid> grid
        = StructuredGridFactory<Grid>::createSimplexGrid
                                        (lower, upper, numElements);
    // To prevent bug in UGGrid destructor
    grids.push_back(grid);

    using LeafGridView = typename Grid::LeafGridView;
    const LeafGridView gridView = grid->leafGridView();

    using FEBasisTest
        = Functions::BernsteinDGRefinedDGBasis<LeafGridView, 1, 3>;
    auto unnormalizedTestSpaces = make_space_tuple<FEBasisTest>(gridView);

    auto unnormalizedInnerProduct =
      make_InnerProduct(unnormalizedTestSpaces,
          make_tuple(
              make_IntegralTerm<0,0,IntegrationType::gradGrad,
                                    DomainOfIntegration::interior>(1., s),
              make_IntegralTerm<0,0,
                                IntegrationType::travelDistanceWeighted,
                                DomainOfIntegration::face>(1., s)));
    double condition;
    if(normalized) {
      auto testSpaces = make_normalized_space_tuple(unnormalizedInnerProduct);
      auto innerProduct
          = replaceTestSpaces(unnormalizedInnerProduct, testSpaces);
      condition = conditionOnFirstElement(innerProduct);
    } else {
      condition = conditionOnFirstElement(unnormalizedInnerProduct);
    }
    std::cout << level << ": " << condition << '\n';
  }
}
