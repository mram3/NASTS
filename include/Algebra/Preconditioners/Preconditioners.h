#pragma once

#include "Matrix/CsrMatrix.h"

#include <vector>

class Preconditioners
{
public:
    //Destructor
    virtual ~Preconditioners() = default;

    virtual void setup
    (
        const CsrMatrix& Mat
    ) 
    = 0;

    virtual void apply
    (
        std::vector<double>& z,
        const std::vector<double>& r
    )
    const = 0;
};