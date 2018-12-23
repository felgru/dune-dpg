// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_FUNCTIONPLOTTER_HH
#define DUNE_DPG_FUNCTIONPLOTTER_HH

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dune/dpg/functions/constraineddiscreteglobalbasisfunction.hh>
#include <dune/dpg/functions/discreteglobalbasisfunction.hh>
#include <dune/dpg/functions/io/vtkrefinedfunctionwriter.hh>
#include <dune/dpg/type_traits.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <algorithm>
#include <string>
#include <type_traits>

namespace Dune {

class FunctionPlotter
{
public:
  FunctionPlotter() = delete;

  FunctionPlotter(std::string filename) : filename(filename) {}

  template<class CoefficientVector, class FEBasis,
      typename std::enable_if<!is_RefinedFiniteElement<FEBasis>::value
                             >::type* = nullptr>
  void plot(const std::string&       functionname,
            const CoefficientVector& coefficientVector,
            const FEBasis&           feBasis,
            unsigned int             subsampling) const
  {
    auto discreteFunction
        = discreteGlobalBasisFunction(feBasis, coefficientVector);
    auto localDiscreteFunction = localFunction(discreteFunction);

    // We need to subsample the grid for plotting, because VTK can only
    // display piecewise polynomials of order 1. Thus we approximate the
    // given discrete function by first order piecewise polynomials
    // living on a finer grid.
    SubsamplingVTKWriter<typename FEBasis::GridView>
        vtkWriter(feBasis.gridView(), Dune::refinementLevels(subsampling));
    vtkWriter.addVertexData(localDiscreteFunction,
                            VTK::FieldInfo(functionname,
                                           VTK::FieldInfo::Type::scalar,
                                           1));
    vtkWriter.write(filename);
  }

  template<class CoefficientVector, class FEBasis,
      typename std::enable_if<is_RefinedFiniteElement<FEBasis>::value
                             >::type* = nullptr>
  void plot(const std::string&       functionname,
            const CoefficientVector& coefficientVector,
            const FEBasis&           feBasis,
            unsigned int             subsampling) const
  {
    // TODO: currently subsampling is ignored for refined finite elements
    std::ofstream file(filename+".vtu");
    VTKRefinedFunctionWriter writer(file);
    writer.writeFunction(feBasis, coefficientVector, functionname);
  }

  template<class CoefficientVector, class FEBasis>
  void plot(const std::string&       functionname,
            const CoefficientVector& x,
            const FEBasis&           feBasis,
            unsigned int             subsampling,
            size_t                   spaceOffset) const
  {
    CoefficientVector coefficientVector(feBasis.size());
    assert(x.size() >= feBasis.size()+spaceOffset);
    std::copy_n(x.begin()+spaceOffset, feBasis.size(),
                coefficientVector.begin());

    plot(functionname, coefficientVector, feBasis, subsampling);
  }

private:
  const std::string filename;
};

} // end namespace Dune

#endif // DUNE_DPG_FUNCTIONPLOTTER_HH
