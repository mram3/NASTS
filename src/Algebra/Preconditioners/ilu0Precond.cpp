#include "Algebra/Preconditioners//IncompleteFamPrecond.h"
#include "Matrix/CsrMatrix.h"

#include <vector>

void ilu0Precond::setup
(
    const CsrMatrix& A
)
{
    matPtr_ = &A;

    const auto& rowPtr = A.pattern()->rowPtr();
    const auto& col = A.pattern()->column();
    const auto& values = A.values();

    luvalues_ = values;

    int n = rowPtr.size() - 1;

    //ILU(0) factorization
    for(int i = 1; i < n; ++i){//we can skip first row because u[0][j] = a[0][j]

        //Loop over off-diagonal elements of row i
        for(int j1 = rowPtr[i]+1; j1 < rowPtr[i+1]; ++j1){
            int k = col[j1];//current column

            //we're looking for lower triangular elements within the row
            //if all lower elements are zero no decomposition needs to be done there
            if(i<=k) break;

            //dividing the lower triangular element
            //with the diagonal element of its own column (row k)
            //Lik = Lik/Dkk
            luvalues_[j1] /= luvalues_[rowPtr[k]];
            
            //Loop over off-diagonal upper diag elements of row k 
            for(int j2 = rowPtr[k]+1; j2<rowPtr[k+1]; ++j2){
                //now going to the row of dividing diagonal and searching for nnz elements
                int target = col[j2];

                //target has to be upper diagonal
                if(target <= k) continue;

                if(target == i){
                    luvalues_[rowPtr[i]] -= luvalues_[j1]*luvalues_[j2];
                }
                else{
                    for(int l = j1 +1; l < rowPtr[i+1]; ++l){
                    //coming back to the original row & checking if there's non zero value in the 
                    //column corresponding to non zero value in row of dividing diagonal
                        if(col[l] == target){
                            luvalues_[l] -= luvalues_[j1]*luvalues_[j2];
                            //if a match is found, we do Aij = Aij - Lik*Ukj
                            break;
                        }
                        else if(col[l] > target){
                            break;
                        }
                    //if a match is not found we do nothing, eliminating the fill in
                    //and preserving the original CSR (sparsity) structure
                    }
                }
            }
        }
    }
}

void ilu0Precond::apply
(
    std::vector<double>& z,
    const std::vector<double>& r
) const
{
    const auto& rowPtr = matPtr_->pattern()->rowPtr();
    const auto& col = matPtr_->pattern()->column();

    int n = r.size();
    
    //ILU(0) Implementation
    //to solve: (LU)z = r
    //Let Uz = y -> solve Ly = r (forward elimination)
    //Then solve U*z = y (Backward substitution)

    z[0] = r[0];
    for(int i = 1; i < n; ++i){//forward elimination
        double sum = 0.0;
        for(int j = rowPtr[i] + 1; j < rowPtr[i+1]; ++j){

            if(i < col[j]) break;

            sum += luvalues_[j]*z[col[j]];

        }
        z[i] = r[i] - sum;
    }

    for(int i = n-1; i >=0; --i){//backward substitution 
        double sum = 0.0;
        for(int j = rowPtr[i] + 1; j < rowPtr[i+1]; ++j){

            if(i > col[j]) continue;

            sum += luvalues_[j] * z[col[j]];
            
        }
        z[i] = (z[i]-sum) / luvalues_[rowPtr[i]];
    }
}