#pragma once

#include "Matrix/FvMatrix.h"
#include "Fields/ScalarField.h"
#include "Fields/VectorField.h"
#include "Type.h"

namespace explicitCalc{

    VectorField grad
    (
        const ScalarField& phi
    );
}