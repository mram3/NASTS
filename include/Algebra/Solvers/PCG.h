#pragma once

#include "Algebra/Solvers/LinearSolvers.h"

class PCG : public LinearSolvers
{
public:
    //inheriting constructor from base class
    using LinearSolvers::LinearSolvers;

    void solve
    (
        const CsrMatrix& A,
        const std::vector<double>& b,
        std::vector<double>& x
    ) const override ;
};