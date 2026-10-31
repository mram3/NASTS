#include "Controls/Time.h"
#include "Config.h"

#include <iostream>

Time::Time
(

)
{
    auto tSettings = Config::access().data()["Time"];
    
    if(tSettings["Scheme"] == "SteadyState"){
        isSteadySolve_ = true;
    }

    if(isSteadySolve_){
        maxIter_ = tSettings["MaxIter"];
        currIter_ = 0;
        dt_ = 1.0;
        currTime_ = 0.0;
    }

    else{
        dt_ = tSettings["dt"];
        start_ = tSettings["StartTime"];
        end_ = tSettings["EndTime"];
        currTime_ = start_;
        currStep_ = 0;
    }

    writeInterval_ = tSettings.value("WriteInterval", 1);
}

bool Time::isRunning
(

)
{
    if(isSteadySolve_){
        if(currIter_ < maxIter_){
            currIter_++;

            std::cout << "Iteration : " << currIter_ << "\n\n";
            return true;
        }

        std::cout << "Max Iterations Reached\n";
        return false;
    }

    if(currTime_ < end_){
        currTime_ += dt_;
        currStep_++;

        std::cout << "Time : " << currTime_ << "\n\n";
        return true;
    }

    std::cout << "End Time reached\n";
    return false;
}

bool Time::write
(

)const
{
    if(writeInterval_ <= 0) return false;

    if(isSteadySolve_){
        return (currIter_ % writeInterval_ == 0);
    }
    else{
        return (currStep_ % writeInterval_ == 0);
    }
}