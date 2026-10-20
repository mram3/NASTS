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
} // namespace implicitAssembly