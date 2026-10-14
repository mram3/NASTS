#include "Matrix/CsrMatrix.h"

#include <vector>

std::vector<double> CsrMatrix::SpMV
(
    const std::vector<double>& x
) const
{
    const auto& rowPtr_ = rowPtr();
    const auto& column_ = column();
    int n = rowPtr().size()-1;

    std::vector<double> result(n, 0.0);

    for(int i = 0; i < n; ++i){
        double sum = 0.0;
        for(int j = rowPtr_[i]; j < rowPtr_[i+1]; ++j){
            sum += values_[j]*x[column_[j]]; 
        }
        result[i] = sum;
    }

    return result;
}