#pragma once

#include "FiniteVolumeMethods/implicitAssembly/DdtSchemes/DdtSchemes.h"

namespace implicitAssembly
{

class SOUEDdt : public DdtSchemes
{
public:

    FvMatrix Ddt
    (
        double rho,
        ScalarField& phi
    ) override;
};

}