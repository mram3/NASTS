#pragma once

#include "Grid/Grid.h"

#include <vector>

class CsrAddress
{
private:
    //CSR address Arrays
    std::vector<int> rowPtr_;
    std::vector<int> column_;

    std::vector<int> oValPtr_;
    std::vector<int> nValPtr_;
public:

    //construct from Grid
    CsrAddress
    (
        const Grid& mesh
    );

    //getters to address

    //for solvers
    const std::vector<int>& rowPtr() const {return rowPtr_; }
    const std::vector<int>& column() const {return column_; }

    //for assembly
    const std::vector<int>& oValPtr() const {return oValPtr_; }
    const std::vector<int>& nValPtr() const {return nValPtr_; }
};