#include "Controls/SteadyControl.h"
#include "Utils/VectorMath.h"
#include "Config.h"

SteadyControl::SteadyControl
(

)
{
    if(Config::access().data()["Time"].contains("SteadyState")){
        reachSteadyState_ = true;
    }
    if(Config::access().data()["Time"]["Scheme"] == "SteadyState"){
        isSteadySolve_ = true;
    }

    if(reachSteadyState_)
    {
        auto sSettings = Config::access().data()["Time"]["SteadyState"];

        ssTol_ = sSettings.value("Tolerance", 1e-6);
    }

    if(isSteadySolve_){
        auto tSettings = Config::access().data()["Time"];

        ssTol_ = tSettings.value("Tolerance", 1e-6);
    }
}

bool SteadyControl::isReached(const ScalarField& phi)
{
    if(!reachSteadyState_ && !isSteadySolve_) return false;

    currRes_ = VectorMath::L2Norm
    (
        phi.scalarVariable,
        phi.oldScalarVariable
    );

    return (currRes_ < ssTol_);
}