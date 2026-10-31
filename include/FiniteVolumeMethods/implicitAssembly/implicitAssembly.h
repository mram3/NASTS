#pragma once

#include "Matrix/FvMatrix.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "FiniteVolumeFields/VectorField.h"
#include "FiniteVolumeFields/SurfaceScalarField.h"

namespace implicitAssembly
{
    FvMatrix laplacian
    (
        double gamma, 
        ScalarField& phi
    );

    FvMatrix Ddt
    (
        double rho,
        ScalarField& phi
    );

    FvMatrix div
    (
        SurfaceScalarField& flux,
        ScalarField& phi
    );
    
} // namespace implicitAssembly