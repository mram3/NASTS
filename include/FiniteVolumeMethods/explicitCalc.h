#pragma once

#include "Matrix/FvMatrix.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "FiniteVolumeFields/VectorField.h"

namespace explicitCalc{

    VectorField grad
    (
        const ScalarField& phi
    );
}