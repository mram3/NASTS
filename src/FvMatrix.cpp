#include "Matrix/FvMatrix.h"

#include <algorithm>

void FvMatrix::addCoeffs
(
    int row_,
    int col_,
    double value_
)
{
    rowsCOO.push_back(row_);
    colCOO.push_back(col_);
    valCOO.push_back(value_);
}

void FvMatrix::addSource
(
    int row,
    double value
)
{
    source[row] += value;
}

struct COO //struct to combine COO arrays to sort them easily
{
    int row;
    int col;
    double val;

    COO(int row_, int col_, double val_) : row(row_), col(col_), val(val_) {}
};


void FvMatrix::compressToCSR
(

)
{
    int n = rowsCOO.size();

    if(n == 0){//empty matrix
        rowPtr.assign(source.size(), 0);
        return;
    }

    std::vector<COO> arrayCOO;
    arrayCOO.reserve(n);

    for(int i = 0; i < rowsCOO.size(); ++i){
        arrayCOO.emplace_back(rowsCOO[i], colCOO[i], valCOO[i]);
    }

    //sorting the COO array
    std::sort(arrayCOO.begin(), arrayCOO.end(), [](COO A, COO B){
        if(A.row == B.row){
            return A.col < B.col;
        }

        else{
            return A.row < B.row;
        }
    });

    //merging and compressing to CSR format
    double sum = 0.0;
    rowPtr.push_back(0);
    for(int i = 0; i < n; ++i){
        if(i != n-1 && arrayCOO[i].row == arrayCOO[i+1].row&& arrayCOO[i].col == arrayCOO[i+1].col)
        {
            sum += arrayCOO[i].val;
        }
        else{
            sum += arrayCOO[i].val;
            col.push_back(arrayCOO[i].col);
            val.push_back(sum);
            sum = 0.0;

            if(i != n-1){
                for(int j = arrayCOO[i].row; j < arrayCOO[i+1].row; ++j){
                    rowPtr.push_back(col.size());
                }
            }
            else{
                for(int j = arrayCOO[i].row; j < source.size(); ++j){
                    rowPtr.push_back(col.size());
                }
            }
        }
    }

    //COO arrays are no longer needed and are wiped from the memory
    rowsCOO.clear();
    colCOO.clear();
    valCOO.clear();
    rowsCOO.shrink_to_fit();
    colCOO.shrink_to_fit();
    valCOO.shrink_to_fit();
}

std::vector<double> FvMatrix::SpMV
(
    const std::vector<double>& x
) const
{
    int n = rowPtr.size() - 1;
    std::vector<double> r(n);

    for(int i = 0; i<n; ++i){
        double sum = 0.0;
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            sum += val[j] * x[col[j]];
        }
        r[i] = sum;
    }

    return r;
}