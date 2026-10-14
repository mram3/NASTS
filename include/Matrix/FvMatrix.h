#pragma once

#include "Matrix/CsrMatrix.h"
#include "Matrix/CsrAddress.h"
#include "Fields/ScalarField.h"

#include <vector>

class FvMatrix{
private:

    //CSR arrays
    ScalarField& scalarField_;
    const CsrAddress& sparsityPattern_;
    CsrMatrix A_;

    std::vector<double> source_; // b term

public:

    //construct from scalar field
    FvMatrix
    (
        ScalarField& sField
    ) :
        scalarField_(sField),
        sparsityPattern_(scalarField_.getMesh().address()),
        A_(sparsityPattern_) //initialize CsrMatrix A
    {
        int numCells = scalarField_.scalarVariable.size();
        source_.assign(numCells, 0.0); //source term initialized
    }

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
    const CsrMatrix& A() const {return A_; }
    const std::vector<double>& source() const {return source_; }
};