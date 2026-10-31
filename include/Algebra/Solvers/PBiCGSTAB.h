#pragma once

#include "Algebra/Solvers/LinearSolvers.h"

class PBiCGSTAB : public LinearSolvers
{
private:
    std::vector<double> r_;//residual vector
    std::vector<double> zr_;//Bi residual vector
    std::vector<double> p_;//search direction
    std::vector<double> zp_;//preconditioned search direction
    std::vector<double> s_, zs_, t_;
    std::vector<double> Ax_; //Ax  = A * x
    std::vector<double> v_;//v = Ap~

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