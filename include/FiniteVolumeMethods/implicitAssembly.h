#pragma once

#include "Matrix/FvMatrix.h"
#include "FiniteVolumeFields/ScalarField.h"

namespace implicitAssembly
{
    FvMatrix laplacian
    (
        double gamma, 
        ScalarField& phi
    );

    FvMatrix EulerDdt
    (
        double rho,
        ScalarField& U
    );

    FvMatrix Ddt
    (
        double rho,
        ScalarField& U
    );
    
} // namespace implicitAssembly