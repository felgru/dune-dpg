// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_ERRORPLOTTER_HH
#define DUNE_DPG_ERRORPLOTTER_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dune/dpg/functions/concepts.hh>
#include <dune/dpg/functions/constraineddiscreteglobalbasisfunction.hh>

#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/istl/bvector.hh>

#include <cmath>
#include <string>
#include <type_traits>

namespace Dune {

class ErrorPlotter
{
  template<class FEBasis, class Vector,
      typename std::enable_if<models<Functions::Concept
                              ::GlobalBasis<typename FEBasis::GridView>,
                            FEBasis>()>::type* = nullptr>
  static auto
  discreteGlobalBasisFunction(const FEBasis& feBasis, const Vector& u) {
    auto uFunction
        = Dune::Functions::makeDiscreteGlobalBasisFunction<double>
              (feBasis, Dune::TypeTree::hybridTreePath(), u);
    return uFunction;
  }

  template<class FEBasis, class Vector,
      typename std::enable_if<models<Functions::Concept::
            ConstrainedGlobalBasis<typename FEBasis::GridView>,
          FEBasis>()>::type* = nullptr>
  static auto
  discreteGlobalBasisFunction(const FEBasis& feBasis, const Vector& u) {
    auto uFunction = Dune::Functions
        ::makeConstrainedDiscreteGlobalBasisFunction<double>(feBasis, u);
    return uFunction;
  }

public:
  ErrorPlotter() = delete;

  ErrorPlotter(std::string filename) : filename(filename) {}

  template<class VectorType, class FEBasis>
  void plot(const std::string& functionname,
            const VectorType&  elementErrors,
            const FEBasis&     feBasis) const
  {
    auto errorFunction = discreteGlobalBasisFunction(feBasis, elementErrors);
    auto localErrorFunction = localFunction(errorFunction);

    VTKWriter<typename FEBasis::GridView>
        vtkWriter(feBasis.gridView(), VTK::nonconforming);
    vtkWriter.addCellData(localErrorFunction,
                          VTK::FieldInfo(functionname,
                                         VTK::FieldInfo::Type::scalar,
                                         1));
    vtkWriter.write(filename);
  }

  template<class EntitySeed, class GridView>
  void plot(const std::string& functionname,
            const std::vector<std::tuple<EntitySeed, double>>& squaredErrors,
            const GridView&    gridView) const
  {
    const auto& grid = gridView.grid();
    const Functions::LagrangeDGBasis<GridView, 0> feBasis(gridView);

    typedef BlockVector<FieldVector<double,1>> VectorType;
    VectorType elementErrors(feBasis.size());

    auto localView = feBasis.localView();
    for(const auto& cellTuple : squaredErrors) {
      const auto e = grid.entity(std::get<0>(cellTuple));
      localView.bind(e);
      elementErrors[localView.index(0)] = std::sqrt(std::get<1>(cellTuple));
    }

    plot(functionname, elementErrors, feBasis);
  }

private:
  const std::string filename;
};

} // end namespace Dune

#endif // DUNE_DPG_ERRORPLOTTER_HH
