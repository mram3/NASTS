#pragma once

#include "Matrix/CsrAddress.h"

#include <vector>

class CsrMatrix{
private:
    const CsrAddress& sparsePattern_;
    std::vector<double> values_;

public:

    //construct from Sparse Pattern
    CsrMatrix(const CsrAddress& pattern) : sparsePattern_(pattern){
        values_.assign(sparsePattern_.column().size(), 0.0);
    }

    //const Overloading 
    //getter methods for preconditioners and solvers (don't ever call non const CsrMatrices there)
    const std::vector<int>& rowPtr() const {return sparsePattern_.rowPtr(); }
    const std::vector<int>& column() const {return sparsePattern_.column(); }
    const std::vector<double>& values() const {return values_; }

    //to fill the matrix entries during matrix assembly
    //FvMatrix probably is the only place a non const CsrMatrix object is required
    std::vector<double>& values() {return values_; }

    //Sparce Matrix Vector Multiplication
    std::vector<double> SpMV
    (
        const std::vector<double>& x
    ) const;
};