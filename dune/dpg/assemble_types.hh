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

enum class LinearIntegrationType {
    valueFunction,           //!< v f where f is an explicit function known at every point of the domain
    gradFunction,            //!< ∇v · f where f is an explicit function known at every point of the domain
    valueValueScattering,    //!< v Scattering where Scattering is $\int_{sphere} collision_kernel density function$
                             //!<              and density function is a FE polynomial
    valueFunctionScattering, //!< v Scattering where Scattering is $\int_{sphere} collision_kernel density function$
                             //!<              and density function has an explicit formula
    normalVectorFunction     //!< v f n · β, where n is the outer unit normal vector and f is an explicit function known at every point of the domain
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
