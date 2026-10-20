#pragma once

#include <vector>

class SparseAddress
{
private:
    //CSR address Arrays
    std::vector<int> rowPtr_;
    std::vector<int> column_;

public:

    //construct from method's mesh using Addressing utility
    SparseAddress
    (
        std::vector<int> row,
        std::vector<int> col
    )
    :   rowPtr_(std::move(row)),
        column_(std::move(col))
    {}

    //getters to address

    //for solvers
    const std::vector<int>& rowPtr() const {return rowPtr_; }
    const std::vector<int>& column() const {return column_; }
};