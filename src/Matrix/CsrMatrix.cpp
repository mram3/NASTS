#include "Matrix/CsrMatrix.h"

#include <vector>

void CsrMatrix::SpMV
(
    const std::vector<double>& x,
    std::vector<double>& Ax
) const
{
    const auto& rowPtr = sparsity_->rowPtr();
    const auto& column = sparsity_->column();
    size_t n = rowPtr.size()-1;

    for(size_t i = 0; i < n; ++i){
        double sum = 0.0;
        for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j){
            sum += values_[j]*x[column[j]]; 
        }
        Ax[i] = sum;
    }
}