#pragma once

#include "FiniteVolumeFields/ScalarField.h"

class SteadyControl
{
private:
    double initRes_;
    double ssTol_;
    double ssRelTol_;
    double currRes_;

    bool reachSteadyState = false;

public:
    //construct from config file
    SteadyControl();

    /// @brief prime convergence check
    /// @return false if converged
    bool isReached(const ScalarField& phi);
};