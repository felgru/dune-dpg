// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RHS_ASSEMBLER_HH
#define DUNE_DPG_RHS_ASSEMBLER_HH

#include <functional>
#include <list>
#include <map>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

#include <boost/mpl/identity.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/transform.hpp>

#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/adapted/array.hpp>
#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/container/set/convert.hpp>
#include <boost/fusion/algorithm/auxiliary/copy.hpp>
#include <boost/fusion/algorithm/transformation/join.hpp>
#include <boost/fusion/algorithm/transformation/transform.hpp>
#include <boost/fusion/algorithm/transformation/zip.hpp>
#include <boost/fusion/algorithm/iteration/accumulate.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/functional/generation/make_fused_procedure.hpp>
#include <boost/fusion/sequence/intrinsic/value_at.hpp>

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
 * \tparam TestSpaces     tuple of test spaces
 */
template<class TestSpaces>
class RhsAssembler
{
public:
  RhsAssembler () = delete;
  /**
   * \brief constructor for RhsAssembler
   *
   * \note For your convenience, use make_RhsAssembler() instead.
   */
  constexpr RhsAssembler (TestSpaces testSpaces)
             : testSpaces(testSpaces)
  { }

  /**
   * \brief Assemble the rhs vector of a DPG system for a given rhs function.
   *
   * Given a tuple of right hand side functions \p volumeTerms,
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
  const TestSpaces& getTestSpaces() const
  { return testSpaces; }

private:
  TestSpaces     testSpaces;
};

/**
 * \brief Creates an RhsAssembler,
 *        deducing the target type from the types of arguments.
 *
 * \param testSpaces     a tuple of test spaces
 */
template<class TestSpaces>
RhsAssembler<TestSpaces>
make_RhsAssembler(TestSpaces testSpaces)
{
  return RhsAssembler<TestSpaces>(testSpaces);
}

template<class TestSpaces>
template <class LinearForm>
void RhsAssembler<TestSpaces>::
assembleRhs(BlockVector<FieldVector<double,1> >& rhs,
            LinearForm& rhsLinearForm)
{
  using namespace boost::fusion;
  using namespace Dune::detail;

  typedef typename std::tuple_element<0,TestSpaces>::type::GridView GridView;
  GridView gridView = std::get<0>(testSpaces).gridView();

  /* set up global offsets */
  size_t globalTestSpaceOffsets[std::tuple_size<TestSpaces>::value];
  const size_t globalTotalTestSize =
      computeOffsets(globalTestSpaceOffsets, testSpaces);

  // set rhs to correct length -- the total number of basis vectors in the basis
  rhs.resize(globalTotalTestSize);

  // Set all entries to zero
  rhs = 0;

  // Views on the FE bases on a single element
  auto testLocalViews     = genericTransformTuple(testSpaces,
                                                  getLocalViewFunctor());

  auto testLocalIndexSets = as_vector(transform(testSpaces,
                                                getLocalIndexSet()));

  for(const auto& e : elements(gridView)) {

    Hybrid::forEach(testLocalViews, applyBind<decltype(e)>(e));

    for_each(zip(testLocalIndexSets, testLocalViews),
             make_fused_procedure(bindLocalIndexSet()));

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> > localRhs;

    rhsLinearForm.bind(testLocalViews);
    rhsLinearForm.getLocalVector(localRhs);

    auto cp = fused_procedure<localToGlobalRHSCopier<decltype(localRhs),
                   typename std::remove_reference<decltype(rhs)>::type> >
                (localToGlobalRHSCopier<decltype(localRhs),
                   typename std::remove_reference<decltype(rhs)>::type>
                                        (localRhs, rhs));

    /* copy every local subvector indexed by a pair of indices from
     * lfIndices exactly once. */
    auto testZip = zip(testLocalViews,
                       testLocalIndexSets,
                       rhsLinearForm.getLocalSpaceOffsets(),
                       globalTestSpaceOffsets);
    typedef
        typename boost::fusion::vector<std::integral_constant<size_t, 0> >
            LFIndices;

    auto lfIndices = LFIndices{};
    for_each(lfIndices,
             localToGlobalRHSCopyHelper<decltype(testZip),
                                        decltype(cp)>
                                       (testZip, cp));
  }
}


template<class TestSpaces>
template <size_t spaceIndex, class ValueType>
void RhsAssembler<TestSpaces>::
applyDirichletBoundary(BlockVector<FieldVector<double,1> >& rhs,
                       const std::vector<bool>& dirichletNodes,
                       const ValueType& value)
{
  static_assert(std::is_arithmetic<ValueType>::value,
                "applyDirichletBoundary not implemented for non arithmetic "
                "boundary data types.");

  const size_t spaceSize =
        std::get<spaceIndex>(testSpaces).size();

  const size_t globalOffset = detail::computeOffset<spaceIndex>(testSpaces);

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
