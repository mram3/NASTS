#pragma once

#include "FiniteVolumeMethods/surfaceInterpolation/surfaceInterpolation.h"

namespace fvinterp
{
class NVF : public surfaceInterpolation
{
public:
    //Destructor
    virtual ~NVF() = default;

    /// @brief Calculates Higher Order High Resolution face value
    /// @param flux - Mass flow across faces
    /// @param phi - Scalar variable being convected
    /// @return - High order face interpolated values
    std::vector<double> interpolateDC
    (
        const SurfaceScalarField& flux,
        const ScalarField& phi
    ) const override;

protected:
    /// @brief Method to calculate limited value defined as pure virtual function here
    /// @param r - ratio of successive gradients 
    /// @return - Limited value
    virtual double limiter
    (
        double r
    ) const
    = 0;
};
}