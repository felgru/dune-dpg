// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_DPG_FUNCTIONS_IO_VTKREFINEDFUNCTIONWRITER_HH
#define DUNE_DPG_FUNCTIONS_IO_VTKREFINEDFUNCTIONWRITER_HH

#include <numeric>
#include <ostream>
#include <string>

#include <dune/dpg/type_traits.hh>
#include <dune/dpg/functions/localindexsetiteration.hh>
#include <dune/grid/io/file/vtk/vtuwriter.hh>

namespace Dune
{
  class VTKRefinedFunctionWriter
  {
    public:
    VTKRefinedFunctionWriter(std::ostream& outstream)
      : writer(outstream, Dune::VTK::ascii, Dune::VTK::unstructuredGrid)
    {};

    template<class FEBasis, class VectorType>
    void writeFunction(const FEBasis& feBasis,
                       const VectorType& data,
                       const std::string& dataName) {
      computeSizes(feBasis);
      writer.beginMain(ncells, npoints);
      writePointData(feBasis, data, dataName);
      writePoints(feBasis);
      writeCells(feBasis);
      writer.endMain();
    }

    private:

    template<class FEBasis>
    void computeSizes(const FEBasis& feBasis) {
      ncells   = 0;
      npoints  = 0;
      ncorners = 0;
      auto localView = feBasis.localView();
      const auto gridView = feBasis.gridView();
      for(const auto& e : elements(gridView)) {
        localView.bind(e);
        const unsigned subcells
            = localView.tree().refinedReferenceElementGridView().size(0);
        ncells += subcells;
        npoints += subcells * 3; // assuming triangles
      }
      // since we have discontinuous data:
      ncorners = npoints;
    }

    template<class LocalView, class VectorType>
    static VectorType
    getLocalDataVector(const LocalView& localView, const VectorType& data)
    {
      VectorType localVector(localView.size());
      copyToLocalVector(data, localVector, localView);
      return localVector;
    }

    template<class FEBasis, class VectorType>
    void writePointData(const FEBasis& feBasis,
                        const VectorType& data,
                        const std::string& dataName)
    {
      writer.beginPointData(dataName);
      {
        std::unique_ptr<VTK::DataArrayWriter> arraywriter
          (writer.makeArrayWriter(dataName, 1, npoints,
                                  VTK::Precision::float32));

        auto localView = feBasis.localView();
        using GridView = typename FEBasis::GridView;
        using ctype = typename GridView::Grid::ctype;
        constexpr static int dim = GridView::dimension;
        const GridView gridView = feBasis.gridView();
        std::vector<FieldVector<double,1>> feValues;
        feValues.reserve(localView.maxSize());
        for(const auto& e : elements(gridView)) {
          localView.bind(e);
          const VectorType localData = getLocalDataVector(localView, data);
          const auto refinementGridView
              = localView.tree().refinedReferenceElementGridView();
          unsigned int subElementOffset = 0;
          localView.resetSubElements();
          for(const auto& eSub : elements(refinementGridView)) {
            localView.bindSubElement(eSub);
            const auto& localFiniteElement = localView.tree().finiteElement();
            const auto referenceElement
                = Dune::referenceElement<ctype, dim>(eSub.type());
            for(int i = 0, iEnd = referenceElement.size(dim); i < iEnd; i++) {
              const auto corner = referenceElement.position(i, dim);
              localFiniteElement.localBasis().evaluateFunction(corner,
                                                               feValues);
              const auto value =
                  std::inner_product(
                    feValues.begin(), feValues.end(),
                    localData.begin() + subElementOffset, 0.);
              arraywriter->write(value);
            }
            if(is_DGRefinedFiniteElement<FEBasis>::value)
              subElementOffset += localFiniteElement.size();
          }
        }
      }
      writer.endPointData();
    }

    template<class FEBasis>
    void writePoints(const FEBasis& feBasis) {
      writer.beginPoints();
      {
        std::unique_ptr<VTK::DataArrayWriter> arraywriter
          (writer.makeArrayWriter("Coordinates", 3, npoints,
                                  VTK::Precision::float32));
        auto localView = feBasis.localView();
        using GridView = typename FEBasis::GridView;
        const GridView gridView = feBasis.gridView();
        for(const auto& e : elements(gridView)) {
          const auto eGeometry = e.geometry();
          localView.bind(e);
          const auto refinementGridView
              = localView.tree().refinedReferenceElementGridView();
          for(const auto& eSub : elements(refinementGridView)) {
            const auto eSubGeometry = eSub.geometry();
            for(int i = 0, iEnd = eSubGeometry.corners(); i < iEnd; i++) {
              const auto corner = eGeometry.global(eSubGeometry.corner(i));
              constexpr int dimw = GridView::dimensionworld;
              int comp = 0;
              for(; comp < dimw; comp++)
                arraywriter->write(corner[comp]);
              for(; comp < 3; comp++)
                arraywriter->write(0);
            }
          }
        }
      }
      writer.endPoints();
    }

    template<class FEBasis>
    void writeCells(const FEBasis& feBasis) {
      using GridView = typename FEBasis::GridView;
      writer.beginCells();
      { // connectivity
        std::unique_ptr<VTK::DataArrayWriter> arraywriter
          (writer.makeArrayWriter("connectivity", 1, ncorners,
                                  VTK::Precision::int32));
        for(unsigned i = 0; i < ncorners; i++)
          arraywriter->write(i);
      }
      { // connectivity
        std::unique_ptr<VTK::DataArrayWriter> arraywriter
          (writer.makeArrayWriter("offsets", 1, ncells,
                                  VTK::Precision::int32));
        int offset = 0;
        auto localView = feBasis.localView();
        const GridView gridView = feBasis.gridView();
        for(const auto& e : elements(gridView)) {
          localView.bind(e);
          const auto refinementGridView
              = localView.tree().refinedReferenceElementGridView();
          for(const auto& eSub : elements(refinementGridView)) {
            const auto eSubGeometry = eSub.geometry();
            offset += eSubGeometry.corners();
            arraywriter->write(offset);
          }
        }
      }
      { // types
        std::unique_ptr<VTK::DataArrayWriter> arraywriter
          (writer.makeArrayWriter("types", 1, ncells, VTK::Precision::uint8));
        auto localView = feBasis.localView();
        const GridView gridView = feBasis.gridView();
        for(const auto& e : elements(gridView)) {
          localView.bind(e);
          const auto refinementGridView
              = localView.tree().refinedReferenceElementGridView();
          for(const auto& eSub : elements(refinementGridView)) {
            int vtktype = VTK::geometryType(eSub.type());
            arraywriter->write(vtktype);
          }
        }
      }
      writer.endCells();
    }

    VTK::VTUWriter writer;
    unsigned ncells;
    unsigned npoints;
    unsigned ncorners;
  };
}

#endif
