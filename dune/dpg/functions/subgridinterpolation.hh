// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_FUNCTIONS_SUBGRIDINTERPOLATION_HH
#define DUNE_DPG_FUNCTIONS_SUBGRIDINTERPOLATION_HH

/** \file
* \brief interpolation between SubGrid and host grid
*/

#include <type_traits>
#include <vector>

#include <dune/dpg/functions/constraineddiscreteglobalbasisfunction.hh>
#include <dune/dpg/functions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

namespace Dune {


template<class SubGridGlobalBasis, class HostGridGlobalBasis,
         class InputVector, class OutputVector>
class SubGridGlobalBasisInterpolator
{
public:
  using Grid = typename SubGridGlobalBasis::GridView::Grid;
  using Element = typename Grid::template Codim<0>::Entity;
  using HostGrid = typename Grid::HostGridType;
  using HostElement = typename HostGrid::template Codim<0>::Entity;

  enum {dim = Grid::dimension};

private:
  struct SubGridFunction {
    using DiscreteGlobalBasisFunction = std::decay_t<decltype(
        discreteGlobalBasisFunction
            (std::declval<SubGridGlobalBasis>(),
             std::declval<InputVector>()))>;

    using HostGeometry = typename HostElement::Geometry;
    using SubGridDomain = typename DiscreteGlobalBasisFunction::Domain;
    using HostDomain = typename HostGeometry::LocalCoordinate;
    using Range = typename DiscreteGlobalBasisFunction::Range;


    struct Traits
    {
       using DomainType = HostDomain;
       using RangeType  = Range;
    };

    SubGridFunction(const SubGridGlobalBasis& subGridBasis,
        const InputVector& subGridCoefficients)
      : subGridFunction(
            discreteGlobalBasisFunction
                (subGridBasis, subGridCoefficients)),
        subGridLocalFunction(localFunction(subGridFunction)),
        hostElement(nullptr)
    {}

    void bind(const Element& element, const HostElement& hostElement) {
      this->hostElement = &hostElement;
      subGridLocalFunction.bind(element);
    }

    void evaluate(const HostDomain& x, Range& y) const {
      const SubGridDomain xsg
        = subGridLocalFunction.localContext().geometry()
          .local(hostElement->geometry().global(x));
      y = subGridLocalFunction(xsg);
    }

    DiscreteGlobalBasisFunction subGridFunction;
    typename DiscreteGlobalBasisFunction::LocalFunction subGridLocalFunction;
    const HostElement* hostElement;
  };


public:
  SubGridGlobalBasisInterpolator(const SubGridGlobalBasis& subGridBasis,
      const InputVector& subGridCoefficients,
      const HostGridGlobalBasis& hostGridBasis,
      OutputVector& hostGridCoefficients)
    : subGridCoefficients(subGridCoefficients),
      hostGridCoefficients(hostGridCoefficients),
      subGridBasis(subGridBasis),
      hostGridBasis(hostGridBasis),
      subGridFunction(subGridBasis, subGridCoefficients)
  {}


  //! Prepare hostgrid vector
  void pre()
  {
    hostGridCoefficients.resize(hostGridBasis.size());
    hostGridCoefficients = 0.0;
  }


  void transfer(const Element& element, const HostElement& hostElement)
  {
    auto subGridLocalView = subGridBasis.localView();
    subGridLocalView.bind(element);
    auto hostGridLocalView = hostGridBasis.localView();
    hostGridLocalView.bind(hostElement);
    const auto& hostGridFiniteElement
        = hostGridLocalView.tree().finiteElement();
    std::vector<double> localHostGridCoefficients;
    subGridFunction.bind(element, hostElement);
    hostGridFiniteElement.localInterpolation().interpolate(subGridFunction,
        localHostGridCoefficients);

    auto hostGridIndexSet = hostGridBasis.localIndexSet();
    hostGridIndexSet.bind(hostGridLocalView);
    for(size_t i = 0, imax = hostGridIndexSet.size();
        i < imax; i++) {
      hostGridCoefficients[hostGridIndexSet.index(i)]
          = localHostGridCoefficients[i];
    }
  }


  //! does nothing
  void post()
  {}


private:
  const InputVector& subGridCoefficients;
  OutputVector& hostGridCoefficients;

  const SubGridGlobalBasis& subGridBasis;
  const HostGridGlobalBasis& hostGridBasis;

  SubGridFunction subGridFunction;
};


/** \brief Interpolate subgrid function to hostgrid
*
* This is a shortcut for creating an (temporary) interpolator
* object and handing it to the transfer method in the subgrid
*/
template<class SubGridGlobalBasis, class HostGridGlobalBasis,
         class InputVector, class OutputVector>
void interpolateFromSubGrid(const SubGridGlobalBasis& subGridBasis,
    const InputVector& subGridFunction,
    const HostGridGlobalBasis& hostGridBasis,
    OutputVector& hostGridFunction)
{
  static_assert(std::is_same<
      typename SubGridGlobalBasis::GridView::Grid::HostGridType,
      typename HostGridGlobalBasis::GridView::Grid>::value,
      "SubGridGlobalBasis and HostGridGlobalBasis are incompatible.");

  using Interpolator
    = SubGridGlobalBasisInterpolator
      <SubGridGlobalBasis, HostGridGlobalBasis, InputVector, OutputVector>;

  Interpolator interpolator(subGridBasis, subGridFunction,
                            hostGridBasis, hostGridFunction);
  subGridBasis.gridView().grid().template transfer<Interpolator>(interpolator);
}

} // namespace Dune

#endif
