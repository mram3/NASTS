#include "Algebra/Preconditioners/GaussFamPrecond.h"
#include "Matrix/CsrMatrix.h"

#include <vector>

void jacobiPrecond::setup
(
    const CsrMatrix& A
)
{
    const auto& rowPtr = A.pattern()->rowPtr();
    const auto& val = A.values();

    int n = rowPtr.size() - 1;

    invDiag_.assign(n, 0.0);

    for(int i = 0; i < n; ++i){
        invDiag_[i] = 1.0/val[rowPtr[i]];
    }
}

void jacobiPrecond::apply
(
    std::vector<double>& z,
    const std::vector<double>& r
) const
{
    int n = r.size();

    for(int i = 0; i < n; ++i){
        z[i] = invDiag_[i] * r[i];
    }
}