#pragma once

#include "FiniteVolumeMethods/surfaceInterpolation/TVD.h"

#include <cmath>

namespace fvinterp
{
class vanLeer : public TVD
{
protected:

double limiter
(
    double r
) const override
{
    return ((r + std::abs(r))/(1 + std::abs(r)));
}
};
}