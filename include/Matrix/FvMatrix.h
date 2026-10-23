#pragma once

#include "Matrix/CsrMatrix.h"
#include "Matrix/SparseAddress.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "nlohmann/json.hpp"

#include <vector>
#include <memory>

class FvMatrix : public CsrMatrix{
private:
    ScalarField& phi_;
    std::vector<double> source_; // b term

    //ptr arrays for Values array in CSR Matrix
    std::vector<int> oValPtr_;//Owner Values Pointer
    std::vector<int> nValPtr_;//Neighbor Values Pointer
public:

    //Destructor 
    virtual ~FvMatrix() = default;

    //construct from the scalar field
    FvMatrix
    (
        ScalarField& phi 
    );

    
    /// construct from scalar field and sparsity address
    FvMatrix(
        ScalarField& phi,
        std::shared_ptr<SparseAddress> address
    );

    /// @brief Diagonal Matrix flag
    bool isDiag = false;

    /// @brief fills face pointers
    void fillValPtrs();

    //the prime matrix assembly methods
    void addDiagCoeffs
    (
        int cellID,
        double value
    );

    void addOffDiagCoeffs
    (
        int faceID,
        double oFlux,
        double nFlux
    );

    void addSource
    (
        int row,
        double value
    );

    //getter method for solvers
    const std::vector<double>& source() const {return source_; }

    //getter method for operators
    std::vector<double>& source(){return source_; }

    //val ptr getters
    const std::vector<int>& oValPtr() const {return oValPtr_; }
    const std::vector<int>& nValPtr() const {return nValPtr_; }
    
    //solve methods
    void solve();
};

//Global Solve methods

void solve(FvMatrix&& eqn);

//Global Operators

FvMatrix& operator== 
(
    FvMatrix& mat, 
    double b
);

FvMatrix&& operator== 
(
    FvMatrix&& mat, 
    double b
);

FvMatrix operator-
(
    FvMatrix&& matA,
    FvMatrix&& matB
);