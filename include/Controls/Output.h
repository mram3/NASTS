#pragma once

#include "FiniteVolumeFields/ScalarField.h"
#include "FiniteVolumeFields/VectorField.h"
#include "Matrix/FvMatrix.h"
#include "Controls/Time.h"

class Output
{
private:
    const FvMesh& fvmesh_;
    const Time& time_;

    std::vector<const ScalarField*> scalarFields_;
    std::vector<const VectorField*> vectorFields_;

public:
    //Constructor
    Output
    (
        const FvMesh& fvmesh,
        const Time& time
    )
    :   fvmesh_(fvmesh),
        time_(time)
    {}

    //to register the fields at create phase
    void registerScalar
    (
        const ScalarField& sF
    )
    {
        scalarFields_.push_back(&sF);
    }

    void registerVector
    (
        const VectorField& vF
    )
    {
        vectorFields_.push_back(&vF);
    }

    //to write the legacy vtk ASCII file 
    void writeVtk
    (
        bool force = false
    ) const;
};