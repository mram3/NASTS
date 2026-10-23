#pragma once

class Time
{
private:
    double dt_;
    double start_;
    double end_;
    double currTime_;

public:

    //Construct from Config file
    Time();

    //Prime time control check
    bool isRunning();

    // getters
    const double currTime() const {return currTime_; }
};