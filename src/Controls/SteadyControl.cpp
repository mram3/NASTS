#include "Controls/SteadyControl.h"
#include "Utils/VectorMath.h"
#include "Config.h"

SteadyControl::SteadyControl
(

)
:   reachSteadyState(Config::access().data()["Time"].contains("SteadyState"))
{
    if(reachSteadyState)
    {
        auto sSettings = Config::access().data()["Time"]["SteadyState"];

        ssTol_ = sSettings.value("Tolerance", 1e-6);
        ssRelTol_ = sSettings.value("RelTol", 0.0);
    }
}

bool SteadyControl::isReached(const ScalarField& phi)
{
    if(!reachSteadyState) return false;

    currRes_ = VectorMath::L2Norm
    (
        phi.scalarVariable,
        phi.oldScalarVariable
    );

    if(currRes_ < ssTol_) return true;

    return false;
}