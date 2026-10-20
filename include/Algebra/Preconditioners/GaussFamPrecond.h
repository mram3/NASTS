#pragma once

#include "Algebra/Preconditioners/Preconditioners.h"

#include <vector>

class jacobiPrecond : public Preconditioners
{
private:
    std::vector<double> invDiag_;

public:

    void setup
    (
        const CsrMatrix& A
    ) override;

    void apply
    (
        std::vector<double>& z, 
        const std::vector<double>& r
    ) const override;
};