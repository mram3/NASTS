#pragma once

#include "Fields/ScalarField.h"

#include <vector>

class FvMatrix{
private:
    // COO arrays
    std::vector<int> rowsCOO;
    std::vector<int> colCOO;
    std::vector<double> valCOO;

    //CSR arrays
    std::vector<int> rowPtr;
    std::vector<int> col;
    std::vector<double> val;

    std::vector<double> source;

public:

    const ScalarField& scalarField;

    //constructor
    FvMatrix(const ScalarField& sField_) : scalarField(sField_)
    {
        int numCells = scalarField.scalarVariable.size();

        //2 coefficients per face * ~5 non zero terms
        rowsCOO.reserve(10 * numCells);
        colCOO.reserve(10 * numCells);
        valCOO.reserve(10 * numCells);

        //row count will be exactly equal to no of cells + 1 to store total non zero terms
        rowPtr.reserve(numCells + 1);

        //2D quad cell in fvm has ~5 non zero terms every row
        col.reserve(5 * numCells);
        val.reserve(5 * numCells);

        source.assign(numCells, 0.0);
    }

    //to dynamically push coefficients to COO arrays
    void addCoeffs
    (
        int row,
        int col,
        double value
    );

    void addSource
    (
        int row,
        double value
    );

    void compressToCSR();

    std::vector<double> SpMV
    (
        const std::vector<double>& x
    ) const;

    const std::vector<double>& getSource() const {return source;}
};