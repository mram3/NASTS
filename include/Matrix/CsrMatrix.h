#pragma once

#include "Matrix/SparseAddress.h"

#include <memory>
#include <vector>

class CsrMatrix{
protected:
    std::shared_ptr<SparseAddress> sparsity_;
    std::vector<double> values_;

public:

    //Destructor
    virtual ~CsrMatrix() = default;

    //construct from Sparse Pattern
    explicit CsrMatrix
    (
        std::shared_ptr<SparseAddress> pattern
    ) 
    :   sparsity_(pattern),
        values_(sparsity_->column().size(), 0.0)
    {}

    //getter methods for preconditioners and solvers (don't ever call non const CsrMatrices there)
    const std::shared_ptr<SparseAddress>& pattern() const {return sparsity_; }
    const std::vector<double>& values() const {return values_; }

    //non const getter for testing
    std::vector<double>& values(){return values_; }
    
    //Sparce Matrix Vector Multiplication
    void SpMV
    (
        const std::vector<double>& x,
        std::vector<double>& Ax
    ) const;
};