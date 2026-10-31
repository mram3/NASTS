#include "FiniteVolumeFields/SurfaceScalarField.h"
#include "Type.h"

#include <cmath>

void SurfaceScalarField::calcFlux
(
    double rho,
    const VectorField& U
)
{
    const auto& faces = U.fvmesh().faces();
    const auto& cells = U.fvmesh().cells();
    int numInFace = U.fvmesh().numInternalfaces();
    int numFaces = faces.size();

    for(int i = 0; i < numInFace; ++i){
        const auto& face = faces[i];

        int oID = face.ownerCellID;
        int nID = face.neighborCellID;

        double w = (face.xOf * face.xOf + face.yOf * face.yOf);
        w /= (face.xON * face.xON + face.yON * face.yON);
        
        w = 1 - std::sqrt(w);

        double u = w * U.vectorVariable[oID].x + (1 - w) * U.vectorVariable[nID].x;
        double v = w * U.vectorVariable[oID].y + (1 - w) * U.vectorVariable[nID].y;

        double Uf = u * face.face_normal_x + v * face.face_normal_y;

        faceValues[i] = rho * Uf * face.face_area;
    }

    for(int b = numInFace; b < numFaces; ++b){
        const auto& face = faces[b];

        int tag = face.physicalTag;

        auto boundCond = U.boundaryConditionMap.at(tag);

        double u = 0.0;
        double v = 0.0;
        if(boundCond.bc == BoundaryConditionType::Dirichlet){
            u = boundCond.value.x;
            v = boundCond.value.y;
        }
        else if(boundCond.bc == BoundaryConditionType::Neumann){
            int oID = face.ownerCellID;

            u = U.vectorVariable[oID].x + boundCond.value.x * face.xOf;
            v = U.vectorVariable[oID].y + boundCond.value.y * face.yOf;
        }

        double Uf = u * face.face_normal_x + v * face.face_normal_y;

        faceValues[b] = rho * Uf * face.face_area;
    }
}