#ifndef DUNE_DPG_ASSEMBLE_TYPES_HH
#define DUNE_DPG_ASSEMBLE_TYPES_HH

#include <utility>

namespace Dune {

//! The domain over which we integrate
enum class DomainOfIntegration {
    interior, //!< integration over the interior of a cell
    face      //!< integration over the faces of a cell
};

//! The type of the integrand in a term of BilinearForm or InnerProduct
enum class IntegrationType {
    valueValue,   //!< v u
    gradValue,    //!< ∇v ·  u
    valueGrad,    //!<  v · ∇u
    gradGrad,     //!< ∇v · ∇u
    normalVector, //!< v u n · β, where n is the outer unit normal vector
    normalSign,
    travelDistanceWeighted
};

//! The type of the integrand in a linear integral term
enum class LinearIntegrationType {
    valueFunction,           //!< v f where f is an explicit function known at every point of the domain
    gradFunction,            //!< ∇v · f where f is an explicit function known at every point of the domain
    normalVectorFunction     //!< v f n · β, where n is the outer unit normal vector and f is an explicit function known at every point of the domain
};

//! How a given function should be evaluated
enum class EvaluationType {
    value,   //!< v
    grad,    //!< ∇v
};

template <class lhsSpaceIndex,
          class rhsSpaceIndex,
          class BilinearTerm>
struct BilinearTermWithIndices
{
  using LhsIndex = lhsSpaceIndex;
  using RhsIndex = rhsSpaceIndex;
  using Term = BilinearTerm;
  Term term;

  BilinearTermWithIndices(Term&& term)
    : term(std::move(term)) {};
};

template <class spaceIndex,
          class LinearTerm>
struct LinearTermWithIndex
{
  using Index = spaceIndex;
  using Term = LinearTerm;
  Term term;

  LinearTermWithIndex(Term&& term)
    : term(std::move(term)) {};
};

} // end namespace Dune

#endif // DUNE_DPG_ASSEMBLE_TYPES_HH
