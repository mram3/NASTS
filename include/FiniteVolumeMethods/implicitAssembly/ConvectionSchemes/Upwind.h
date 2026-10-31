#pragma once

#include "FiniteVolumeMethods/implicitAssembly/ConvectionSchemes/ConvectionSchemes.h"

namespace implicitAssembly
{
class Upwind : public ConvectionScheme
{
public: 
    FvMatrix div
    (
        SurfaceScalarField& flux,
        ScalarField& phi
    ) const override;
};

}