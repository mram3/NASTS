#include "FiniteVolumeMethods/implicitAssembly/ConvectionSchemes/Upwind.h"

namespace implicitAssembly
{

FvMatrix Upwind::div
(
    SurfaceScalarField& flux,
    ScalarField& phi
) const
{
    //initializing full spatial matrix array
    FvMatrix mat(phi);

    const auto& cells = phi.fvmesh().cells();
    const auto& faces = phi.fvmesh().faces();
    int numInFaces = phi.fvmesh().numInternalfaces();
    int n = faces.size();

    // brief Upwind Differencing Scheme
    for(int i = 0; i < numInFaces; ++i){//looping through internal faces
        const auto& face = faces[i];

        int oID = face.ownerCellID;
        int nID = face.neighborCellID;

        double M = flux.faceValues[i];
        if(M > 0){
            mat.addDiagCoeffs(oID, M);
            mat.addOffDiagCoeffs(i, 0.0, -M);
        }
        else{
            mat.addDiagCoeffs(nID, -M);
            mat.addOffDiagCoeffs(i, M, 0.0);
        }
    }

    for(int b = numInFaces; b < n; ++b){
        const auto& face = faces[b];

        int oID = face.ownerCellID;
        int tag = face.physicalTag;

        auto id = phi.boundaryConditionMap.find(tag);

        auto bcType = id->second.bc;
        auto bcValue = id->second.value;

        double M = flux.faceValues[b];
        if(bcType == BoundaryConditionType::Dirichlet){

            if(M < 0){
                double s = M * bcValue;
                mat.addSource(oID, -s);
            }
            else{
                mat.addDiagCoeffs(oID, M);
            }
        }

        else if(bcType == BoundaryConditionType::Neumann){

            if(M > 0){
                mat.addDiagCoeffs(oID, M);
            }
        }
    }

    return mat;
}
};