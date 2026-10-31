#pragma once

#include "FiniteVolumeFields/ScalarField.h"

class SteadyControl
{
private:
    double ssTol_;
    double currRes_;

    bool reachSteadyState_ = false;
    bool isSteadySolve_ = false;

public:
    //construct from config file
    SteadyControl();

    /// @brief prime convergence check
    /// @return false if converged
    bool isReached(const ScalarField& phi);
};