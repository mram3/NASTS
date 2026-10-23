#include "Controls/Time.h"
#include "Config.h"

#include <iostream>

Time::Time
(

)
{
    auto tSettings = Config::access().data()["Time"];
    
    dt_ = tSettings["dt"];
    start_ = tSettings["StartTime"];
    end_ = tSettings["EndTime"];
    currTime_ = start_;
}

bool Time::isRunning
(

)
{
    if(currTime_ < end_){
        currTime_ += dt_;
        return true;
    }

    std::cout << "End Time reached\n";
    return false;
}