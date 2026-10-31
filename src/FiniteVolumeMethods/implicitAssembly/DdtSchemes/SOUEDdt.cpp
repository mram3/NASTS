#include "FiniteVolumeMethods/implicitAssembly/DdtSchemes/SOUEDdt.h"
#include "FiniteVolumeMethods/implicitAssembly/implicitAssembly.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "Matrix/FvMatrix.h"
#include "Config.h"

#include <iostream>

namespace implicitAssembly
{

FvMatrix SOUEDdt::Ddt
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

        double rhoVbyDt = rho * cell.cell_volume/dt;

        double ap = 1.5 * rhoVbyDt;
        double b = 2 * rhoVbyDt * phi.oldScalarVariable[cell.cellID]
                 - 0.5 * rhoVbyDt * phi.tooOldScalarVariable[cell.cellID];

        Mat.addDiagCoeffs(cell.cellID, ap);

        Mat.addSource(cell.cellID, b);
    }

    return Mat;
}
}