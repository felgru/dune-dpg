#ifndef DUNE_DPG_ASSEMBLE_TYPES_HH
#define DUNE_DPG_ASSEMBLE_TYPES_HH

namespace Dune {

enum class DomainOfIntegration {
    interior,
    face,
    subface
};

enum class IntegrationType {
    valueValue,
    gradValue,
    valueGrad,
    gradGrad
};

enum class SpaceType {
    test,
    solution
};

} // end namespace Dune

#endif // DUNE_DPG_ASSEMBLE_TYPES_HH
