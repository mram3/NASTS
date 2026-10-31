#pragma once

#include "FieldSpace.h"
#include "FiniteVolumeFields/VectorField.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "Type.h"
#include "nlohmann/json.hpp"

#include <vector>

class SurfaceScalarField : public FieldSpace
{
public:
    //Storage for face values like flux through faces etc.,
    std::vector<double> faceValues;

    SurfaceScalarField
    (
        std::string& fieldName,
        const FvMesh& fvmesh
    )
    :   FieldSpace(fieldName, fvmesh),
        faceValues(fvmesh.faces().size(), 0.0)     
    {}

    void calcFlux
    (
        double rho,
        const VectorField& U
    );
};