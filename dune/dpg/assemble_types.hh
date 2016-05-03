#ifndef DUNE_DPG_ASSEMBLE_TYPES_HH
#define DUNE_DPG_ASSEMBLE_TYPES_HH

namespace Dune {

enum class DomainOfIntegration {
    interior, //!< integration over the interior of a cell
    face      //!< integration over the faces of a cell
};

enum class IntegrationType {
    valueValue,   //!< v u
    gradValue,    //!< ∇v ·  u
    valueGrad,    //!<  v · ∇u
    gradGrad,     //!< ∇v · ∇u
    normalVector, //!< v u n · β, where n is the outer unit normal vector
    normalSign,
    travelDistanceWeighted
};

enum class EvaluationType {
    value,   //!< v
    grad,    //!< ∇v
};

enum class SpaceType {
    test,     //!< a test space
    solution  //!< a solution space
};

//! The saddlepoint formulation of a DPG system.
struct SaddlepointFormulation {};
//! The optimal test function formulation of a DPG system.
struct DPGFormulation {};

} // end namespace Dune

#endif // DUNE_DPG_ASSEMBLE_TYPES_HH
