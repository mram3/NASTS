#include "FiniteVolumeMethods/implicitAssembly/DdtSchemes/EulerDdt.h"
#include "FiniteVolumeMethods/implicitAssembly/implicitAssembly.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "Matrix/FvMatrix.h"
#include "Config.h"

namespace implicitAssembly
{

FvMatrix EulerDdt::Ddt
(
    double rho,
    ScalarField& phi
)
{
    FvMatrix Mat(phi, phi.fvmesh().fvDiagAddress());
    Mat.isDiag = true;

    const auto& cells = phi.fvmesh().cells();

    double dt = Config::access().data()["Time"]["dt"];
    for(const auto& cell : cells){

        double ap = rho*cell.cell_volume/dt;
        double b = ap * phi.oldScalarVariable[cell.cellID];

        Mat.addDiagCoeffs(cell.cellID, ap);

        Mat.addSource(cell.cellID, b);
    }

    return Mat;
}

} // namespace implicitAssembly