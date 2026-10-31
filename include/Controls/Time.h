#pragma once

class Time
{
private:
    //for unsteady solving
    double dt_;
    double start_;
    double end_;
    double currTime_;

    //for steady state solving
    bool isSteadySolve_ = false;
    int maxIter_;
    int currIter_;

    // I/O control variables
    int writeInterval_;
    int currStep_;

public:

    //Construct from Config file
    Time();

    //Prime time control check
    bool isRunning();

    // getters
    const double currTime() const {return currTime_; }
    const int currIter() const {return currIter_; }
    const int currStep() const {return currStep_; }

    //Flag for solving steady state
    bool isSteadySolve() const {return isSteadySolve_; }

    //Output control method
    bool write() const;
};