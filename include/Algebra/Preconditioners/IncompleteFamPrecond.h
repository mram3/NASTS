#pragma once

#include "Algebra/Preconditioners/Preconditioners.h"

#include <vector>

class ilu0Precond : public Preconditioners
{
public:
    void setup
    (
        const CsrMatrix& A
    ) override;

    void apply
    (
        std::vector<double>& z, 
        const std::vector<double>& r
    )const override;

private:
    std::vector<double> luvalues_;
    const CsrMatrix* matPtr_ = nullptr;
};

class ic0Precond : public Preconditioners
{
public:
    void setup
    (
        const CsrMatrix& A
    ) override;

    void apply
    (
        std::vector<double>& z, 
        const std::vector<double>& r
    )const override;

private:
    std::vector<double> lltvalues_;
    const CsrMatrix* matPtr_ = nullptr;
};