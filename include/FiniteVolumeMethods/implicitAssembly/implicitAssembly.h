#pragma once

#include "Matrix/FvMatrix.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "FiniteVolumeFields/VectorField.h"

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
        double rho,
        ScalarField& phi,
        VectorField& U
    );
    
} // namespace implicitAssembly