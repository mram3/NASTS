#include "FiniteVolumeMethods/implicitAssembly/DdtSchemes/SteadyStateDdt.h"
#include "FiniteVolumeMethods/implicitAssembly/implicitAssembly.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "Matrix/FvMatrix.h"
#include "Config.h"

namespace implicitAssembly
{

FvMatrix SteadyStateDdt::Ddt
(
    double rho,
    ScalarField& phi
)
{
    FvMatrix Mat(phi, phi.fvmesh().fvDiagAddress());
    Mat.isDiag = true;

    return Mat;
}
}