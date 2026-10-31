#pragma once

#include "FiniteVolumeMethods/implicitAssembly/DdtSchemes/DdtSchemes.h"

namespace implicitAssembly
{

class EulerDdt : public DdtSchemes
{
public:

    FvMatrix Ddt
    (
        double rho,
        ScalarField& phi
    ) override;
};

}