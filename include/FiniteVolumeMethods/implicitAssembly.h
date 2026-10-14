#pragma once

#include "Matrix/FvMatrix.h"
#include "Fields/ScalarField.h"
#include "Grid/Grid.h"

namespace implicitAssembly
{
    FvMatrix laplacian
    (
        double gamma, 
        ScalarField& phi
    );
} // namespace implicitAssembly
