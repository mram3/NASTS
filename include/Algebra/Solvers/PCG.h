#pragma once

#include "Algebra/Solvers/LinearSolvers.h"

#include <vector>

class PCG : public LinearSolvers
{
private:
    std::vector<double> r_;//residual vector
    std::vector<double> z_;//preconditioned residual
    std::vector<double> p_;//search direction
    std::vector<double> Ax_;// A * x
    std::vector<double> Ap_;// A * p

public:
    //inheriting constructor from base class
    using LinearSolvers::LinearSolvers;

    SolverPerf solve
    (
        const CsrMatrix& A,
        const std::vector<double>& b,
        std::vector<double>& x
    ) override ;
};