// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_FUNCTIONPLOTTER_HH
#define DUNE_DPG_FUNCTIONPLOTTER_HH

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dune/dpg/functions/concepts.hh>
#include <dune/dpg/functions/constraineddiscreteglobalbasisfunction.hh>

#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <type_traits>

namespace Dune {

class FunctionPlotter
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
  FunctionPlotter() = delete;

  FunctionPlotter(std::string filename) : filename(filename) {}

  template<class VectorType, class FEBasis>
  void plot(const std::string& functionname,
            const VectorType&  u,
            const FEBasis&     feBasis,
            unsigned int       subsampling) const
  {
    //////////////////////////////////////////////////////////////////////////
    //  Make a discrete function from the FE basis and the coefficient vector
    //////////////////////////////////////////////////////////////////////////

    auto uFunction = discreteGlobalBasisFunction(feBasis, u);
    auto localUFunction = localFunction(uFunction);

    //////////////////////////////////////////////////////////////////////////
    //  Write result to VTK file
    //  We need to subsample, because VTK cannot natively display real
    //  second-order functions
    //////////////////////////////////////////////////////////////////////////

    SubsamplingVTKWriter<typename FEBasis::GridView>
        vtkWriter(feBasis.gridView(), Dune::refinementLevels(subsampling));
    vtkWriter.addVertexData(localUFunction,
                            VTK::FieldInfo(functionname,
                                           VTK::FieldInfo::Type::scalar,
                                           1));
    vtkWriter.write(filename);
  }

  template<class VectorType, class FEBasis>
  void plot(const std::string& functionname,
            const VectorType&  x,
            const FEBasis&     feBasis,
            unsigned int       subsampling,
            size_t             spaceOffset) const
  {
    VectorType u(feBasis.size());
    assert(x.size() >= feBasis.size()+spaceOffset);
    std::copy(x.begin()+spaceOffset, x.begin()+spaceOffset+feBasis.size(),
              u.begin());

    plot(functionname, u, feBasis, subsampling);
  }

private:
  const std::string filename;
};

} // end namespace Dune

#endif // DUNE_DPG_FUNCTIONPLOTTER_HH
