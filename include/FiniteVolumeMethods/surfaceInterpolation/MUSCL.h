#pragma once

#include "FiniteVolumeMethods/surfaceInterpolation/TVD.h"

#include <cmath>
#include <algorithm>

namespace fvinterp
{
class MUSCL : public TVD
{
protected:

double limiter
(
    double r
)const override
{
    return std::max(0.0, std::min(2.0 * r, std::min(0.5 * r + 0.5, 2.0)));
}
};
}