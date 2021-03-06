// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_ERRORPLOTTER_HH
#define DUNE_DPG_ERRORPLOTTER_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dune/dpg/functions/discreteglobalbasisfunction.hh>

#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/istl/bvector.hh>

#include <cmath>
#include <string>
#include <type_traits>

namespace Dune {

class ErrorPlotter
{
  using CoefficientVector = BlockVector<FieldVector<double,1>>;

public:
  ErrorPlotter() = delete;

  ErrorPlotter(std::string filename) : filename(filename) {}

  template<class EntitySeed, class GridView>
  void plot(const std::string& functionname,
            const std::vector<std::tuple<EntitySeed, double>>& squaredErrors,
            const GridView&    gridView) const
  {
    const auto& grid = gridView.grid();
    const Functions::LagrangeDGBasis<GridView, 0> feBasis(gridView);

    CoefficientVector elementErrors(feBasis.size());

    auto localView = feBasis.localView();
    for(const auto& [entitySeed, squaredError] : squaredErrors) {
      const auto e = grid.entity(entitySeed);
      localView.bind(e);
      elementErrors[localView.index(0)] = std::sqrt(squaredError);
    }

    plot(functionname, elementErrors, feBasis);
  }

private:
  template<class FEBasis>
  void plot(const std::string& functionname,
            const CoefficientVector& elementErrors,
            const FEBasis& feBasis) const
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

  const std::string filename;
};

} // end namespace Dune

#endif // DUNE_DPG_ERRORPLOTTER_HH
