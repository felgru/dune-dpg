// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RHS_ASSEMBLER_HH
#define DUNE_DPG_RHS_ASSEMBLER_HH

#include <array>
#include <list>
#include <map>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

#include <boost/hana.hpp>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include "assemble_helper.hh"
#include "assemble_types.hh"
#include "bilinearform.hh"
#include "innerproduct.hh"
#include "quadrature.hh"
#include "linearform.hh"
#include "localevaluation.hh"

namespace Dune {

/**
 * \brief This constructs the right hand side vector of a DPG system.
 *
 * \tparam TestSpacesPtr  shared_ptr to a tuple of test spaces
 */
template<class TestSpacesPtr>
class RhsAssembler
{
public:
  using TestSpaces = typename TestSpacesPtr::element_type;

  RhsAssembler () = delete;
  /**
   * \brief constructor for RhsAssembler
   *
   * \note For your convenience, use make_RhsAssembler() instead.
   */
  constexpr RhsAssembler (const TestSpacesPtr& testSpaces)
             : testSpaces(testSpaces)
  { }

  /**
   * \brief Assemble the rhs vector of a DPG system for a given rhs function.
   *
   * Given a LinearForm \p rhsLinearForm,
   * this assembles the right hand side vector of the corresponding
   * DPG system. \p rhs will be overwritten by this function.
   *
   * \param[out] rhs           the rhs vector of the DPG system
   * \param[in]  rhsLinearForm the linear form describing the rhs
   */
  template <class LinearForm>
  void assembleRhs(BlockVector<FieldVector<double,1> >& rhs,
                   LinearForm& rhsLinearForm);

  /**
   * \brief Apply Dirichlet boundary values to the rhs vector
   *
   * \param[in,out] rhs         the rhs vector of the DPG system
   * \param[in] dirichletNodes  true marks the dofs in the Dirichlet boundary
   * \param[in] value           the Dirichlet boundary value
   * \tparam spaceIndex  the index of the test space on which we apply
   *                     the boundary data
   * \tparam ValueType   we take either constants or functions for \p value
   */
  template <size_t spaceIndex, class ValueType>
  void applyDirichletBoundary(BlockVector<FieldVector<double,1> >& rhs,
                              const std::vector<bool>& dirichletNodes,
                              const ValueType& value);

  /**
   * \brief Does exactly what it says on the tin.
   */
  TestSpacesPtr getTestSpaces() const
  { return testSpaces; }

private:
  TestSpacesPtr testSpaces;
};

/**
 * \brief Creates an RhsAssembler,
 *        deducing the target type from the types of arguments.
 *
 * \param testSpaces     a shared_ptr to a tuple of test spaces
 */
template<class TestSpacesPtr>
RhsAssembler<TestSpacesPtr>
make_RhsAssembler(const TestSpacesPtr& testSpaces)
{
  return RhsAssembler<TestSpacesPtr>(testSpaces);
}

template<class TestSpacesPtr>
template <class LinearForm>
void RhsAssembler<TestSpacesPtr>::
assembleRhs(BlockVector<FieldVector<double,1> >& rhs,
            LinearForm& rhsLinearForm)
{
  using namespace Dune::detail;

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(*testSpaces).gridView();

  /* set up global offsets */
  std::array<size_t,std::tuple_size<TestSpaces>::value> globalTestSpaceOffsets;
  const size_t globalTotalTestSize =
      computeOffsets(globalTestSpaceOffsets, *testSpaces);

  // set rhs to correct length -- the total number of basis vectors in the basis
  rhs.resize(globalTotalTestSize);

  // Set all entries to zero
  rhs = 0;

  // Views on the FE bases on a single element
  auto testLocalViews     = getLocalViews(*testSpaces);

  for(const auto& e : elements(gridView)) {

    bindLocalViews(testLocalViews, e);

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> > localRhs;

    rhsLinearForm.bind(testLocalViews);
    rhsLinearForm.getLocalVector(localRhs);

    // TODO: We might copy zero blocks that could be avoided.
    namespace hana = boost::hana;
    auto lfIndices = hana::to<hana::tuple_tag>(
        hana::make_range(hana::int_c<0>,
              hana::int_c<std::tuple_size<TestSpaces>::value>));

    /* copy every local subvector indexed by a pair of indices from
     * LFIndices exactly once. */
    copyLocalToGlobalVector<decltype(lfIndices)>(
        localRhs,
        rhs,
        testLocalViews,
        rhsLinearForm.getLocalSpaceOffsets(),
        globalTestSpaceOffsets);
  }
}


template<class TestSpacesPtr>
template <size_t spaceIndex, class ValueType>
void RhsAssembler<TestSpacesPtr>::
applyDirichletBoundary(BlockVector<FieldVector<double,1> >& rhs,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& value)
{
  static_assert(std::is_arithmetic<ValueType>::value,
                "applyDirichletBoundary not implemented for non arithmetic "
                "boundary data types.");

  const size_t spaceSize =
        std::get<spaceIndex>(*testSpaces).size();

  const size_t globalOffset = detail::computeOffset<spaceIndex>(*testSpaces);

  // Set Dirichlet values
  for (size_t i=0; i<spaceSize; i++)
  {
    if (dirichletNodes[i])
    {
      /* TODO: Needs adaptation when value is a function. */
      rhs[globalOffset+i] = value;
    }
  }
}

} // end namespace Dune

#endif // DUNE_DPG_RHS_ASSEMBLER_HH
