#include "FiniteVolumeMethods/surfaceInterpolation/TVD.h"
#include "FiniteVolumeFields/VectorField.h"
#include "FiniteVolumeMethods/explicitCalc/explicitCalc.h"

#include <vector>

namespace fvinterp
{
/// @brief 12. High Resolution Scheme
///        The FVM in CFD by F. Moukalled, L. Mangani, M. Darwish
///        12.4 The TVD Framework
std::vector<double> TVD::interpolateDC
(
    const SurfaceScalarField& flux,
    const ScalarField& phi
) const
{
    //brief Implementation of Higher Order Schemes in Unstructured Grid Systems
    //      Ref:- 12.6 HR Schemes in Unstructured Grid Systems - Deferred Correction
    //      The FVM in CFD by F. Moukalled, L. Mangani, M. Darwish

    VectorField gradPhi = explicitCalc::grad(phi);

    const auto& faces = phi.fvmesh().faces();
    const auto& cells = phi.fvmesh().cells();
    const auto numInface = phi.fvmesh().numInternalfaces();
    int n = faces.size();

    std::vector<double> phiHO(n, 0.0);
    double phiD = 0.0;
    double phiC = 0.0;
    double gradCdotD = 0.0;

    //To find rf = (phiC - phiU)/(phiD - phiC)
    //where the virtual far upwind node is found by phiU = phiD - 2 * grad(phiC) dot d(CD)\
    //by some rearrangement rf = [2*(grad(phiC)dot d(CD))/(phiD - phiC)] - 1
    for(int i = 0; i < numInface; ++i){
        const auto& face = faces[i];

        if(flux.faceValues[i] > 0){
            phiC = phi.scalarVariable[face.ownerCellID];
            phiD = phi.scalarVariable[face.neighborCellID];

            gradCdotD = gradPhi.vectorVariable[face.ownerCellID].x * face.xON
                        + gradPhi.vectorVariable[face.ownerCellID].y * face.yON;
            
        }
        else{
            phiD = phi.scalarVariable[face.ownerCellID];
            phiC = phi.scalarVariable[face.neighborCellID];

            gradCdotD = - gradPhi.vectorVariable[face.neighborCellID].x * face.xON
                        - gradPhi.vectorVariable[face.neighborCellID].y * face.yON;
            
        }

        double denom = phiD - phiC;
        double r = 0.0;
        double psir = 0.0;

        if(std::abs(gradCdotD) >= 1000 * std::abs(denom)){//to safeguard against denom = 0
            double signNum = gradCdotD > 0 ? 1 : -1;
            double signDen = denom > 0 ? 1 : -1;

            r = 2 * 1000 * signNum * signDen - 1;
        }
        else{
            r = 2 * (gradCdotD/denom) - 1;
        }
        
        psir = limiter(r);
        
        phiHO[i] = phiC + 0.5 * psir * (phiD - phiC);
    }

    return phiHO;
}
}