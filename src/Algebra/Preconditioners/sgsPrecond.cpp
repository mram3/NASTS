#include "Algebra/Preconditioners/GaussFamPrecond.h"
#include "Matrix/CsrMatrix.h"

#include <vector>

void sgsPrecond::setup
(
    const CsrMatrix& A
)
{
    //storing a ref to matrix for apply method
    matPtr_ = &A;

    const auto& rowPtr = A.pattern()->rowPtr();
    const auto& val = A.values();

    int n = rowPtr.size() - 1;

    invDiag_.assign(n, 0.0);

    for(int i = 0; i < n; ++i){
        invDiag_[i] = 1.0/val[rowPtr[i]];
    }
}

void sgsPrecond::apply
(
    std::vector<double>& z,
    const std::vector<double>& r
) const
{
    const auto& rowPtr = matPtr_->pattern()->rowPtr();
    const auto& col = matPtr_->pattern()->column();
    const auto& values = matPtr_->values();

    int n = r.size();

    //Symmetric Gauss Seidel implementation
    //(D+L)*D^-1*(D+U) z = r
    //Let D^-1 * (D+U) z = y -> (D+L) * y = r -> forward substitution
    //then (D+U) z = D * y -> backward elimination to find z

    //forward substitution 
    for(int i = 0; i < n; ++i){
        double sum = 0.0;
        for(int j = rowPtr[i] + 1; j<rowPtr[i+1]; ++j){
            if(i > col[j]){
                sum += values[j] * z[col[j]];
            }
        }
        z[i] = (r[i] - sum) * invDiag_[i];
    }

    //backward elimination
    for(int i = n-1; i >=0; --i){
        double sum = 0.0;
        for(int j = rowPtr[i] + 1; j < rowPtr[i+1]; ++j){
            if(i < col[j]){
                sum += values[j] * z[col[j]];
            }
        }
        z[i] = z[i]- (sum * invDiag_[i]);
    }
}