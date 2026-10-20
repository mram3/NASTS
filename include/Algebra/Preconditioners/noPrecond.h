#pragma once

#include "Algebra/Preconditioners/Preconditioners.h"

class noPrecond : public Preconditioners
{
public:
    void setup
    (
        const CsrMatrix& A
    ) override
    {}

    void apply
    (
        std::vector<double>& z,
        const std::vector<double>& r
    ) const override
    {
        z = r;
    }
};