#pragma once

#include "FiniteVolumeFields/ScalarField.h"
#include "FiniteVolumeFields/VectorField.h"
#include "FiniteVolumeFields/SurfaceScalarField.h"

#include <vector>

namespace fvinterp
{
class surfaceInterpolation
{
public:
    //Destructor
    virtual ~surfaceInterpolation() = default;

    /// @brief Calculates higher order face value for Deferred correction method
    /// @param flux - Flux through faces
    /// @param phi - Scalar Field
    /// @return Field values at the face
    virtual std::vector<double> interpolateDC
    (
        const SurfaceScalarField& flux,
        const ScalarField& phi
    ) const
    = 0;
};
}