#pragma once

#include "FiniteVolumeMethods/implicitAssembly/ConvectionSchemes/ConvectionSchemes.h"
#include "FiniteVolumeMethods/implicitAssembly/ConvectionSchemes/Upwind.h"
#include "FiniteVolumeMethods/surfaceInterpolation/surfaceInterpolation.h"

#include <memory>

namespace implicitAssembly
{
class DeferredCorrection : public ConvectionScheme
{
private:
    Upwind upwind_;
    std::unique_ptr<fvinterp::surfaceInterpolation> surfInterp_;

public:

    /// @brief Constructor
    /// @param surfInterp - Pointer to specific HO face interpolation scheme
    DeferredCorrection
    (
        std::unique_ptr<fvinterp::surfaceInterpolation> surfInterp
    )
    : surfInterp_
    (
        std::move(surfInterp)
    )
    {}

    /// @brief - Arranges Deferred correction Matrix Eqn for HO schemes
    /// @param flux - Flux through face
    /// @param phi - Scalar transported variable
    /// @return - Assembled fvMatrix to main app. loop
    FvMatrix div
    (
        SurfaceScalarField& flux,
        ScalarField& phi
    ) const override;
};

}