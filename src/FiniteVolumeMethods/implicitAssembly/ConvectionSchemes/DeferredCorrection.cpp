#include "FiniteVolumeMethods/implicitAssembly/ConvectionSchemes/DeferredCorrection.h"

#include <vector>
#include <iostream>

namespace implicitAssembly
{

FvMatrix DeferredCorrection::div
(
    SurfaceScalarField& flux,
    ScalarField& phi
) const
{
    // Assemble the Low Order Matrix first
    FvMatrix mat = upwind_.div(flux, phi);

    //High order face interpolation based on previous iterations
    std::vector<double> phiHO = surfInterp_->interpolateDC(flux, phi);

    const auto& faces = phi.fvmesh().faces();
    const auto& cells = phi.fvmesh().cells();
    int numInFace = phi.fvmesh().numInternalfaces();

    for(int i = 0; i < numInFace; ++i){//adding deferred correction to the source term
        const auto& face = faces[i];

        int oID = face.ownerCellID;
        int nID = face.neighborCellID;

        //if(std::isnan(phiHO[i])||std::isinf(phiHO[i])){
        //    std::cerr << "phi HO is Nan/Inf at face " << i << "\n";
        //    std::cerr << "Owner : " << oID << " | Neighbor : " << nID << "\n";
        //}
        //calculate upwind value based on previous phi value
        double phiLO = flux.faceValues[i] > 0 ? phi.scalarVariable[oID] :
                                                    phi.scalarVariable[nID];

        double defCorr = - flux.faceValues[i] * (phiHO[i] - phiLO);

        //if(std::isnan(defCorr) || std::isinf(defCorr)){
        //    std::cerr << "defCorr is NaN/Inf at face : " << i << "\n";
        //    std::cerr << "Flux : " << flux.faceValues[i] 
        //              << " | phiHO : " << phiHO[i] 
        //              << " | phiLO : " << phiLO << "\n";
        //}
        mat.addSource(oID, defCorr);
        mat.addSource(nID, -defCorr);
    }

    return mat;
}
}