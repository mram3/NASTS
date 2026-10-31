#include "Algebra/Preconditioners//IncompleteFamPrecond.h"
#include "Matrix/CsrMatrix.h"

#include <vector>
#include <cmath>

void ic0Precond::setup
(
    const CsrMatrix& A
)
{
    matPtr_ = &A;

    const auto& rowPtr = A.pattern()->rowPtr();
    const auto& col = A.pattern()->column();
    const auto& values = A.values();

    lltvalues_ = values;

    int n = rowPtr.size() - 1;

    //IC(0) factorization
    for(int i = 0; i < n; ++i){

        //to store the dot within the row i before diag col
        double sumD = lltvalues_[rowPtr[i]];//Lii
        
        for(int j1 = rowPtr[i]+1; j1 < rowPtr[i+1]; ++j1){
            int k = col[j1];

            //looking only at lower diagonals
            if(k > i) break;

            //to store dot prod between row i and row k before col k
            double sumL = lltvalues_[j1];//Lik
            for(int j2 = rowPtr[k]+1; j2 < rowPtr[k+1]; ++j2){

                int target = col[j2];//nnz val at row k

                //if target col is after col k - break
                if(target > k) break;

                //coming back to row i to search for corresponding nnz at col target
                for(int j3 = rowPtr[i] + 1; j3 < j1; ++j3){
                    if(col[j3] == target){
                        //if there's a match we do dot prod
                        sumL -= lltvalues_[j2]*lltvalues_[j3];//Lktarget * Litarget
                        break;
                    }
                    else if(col[j3] > target){
                        //early break if there's no nnz value before col target
                        break;
                    }

                //if no match is found we do nothing (eliminating the fill in)
                }
            }

            lltvalues_[j1] = sumL/lltvalues_[rowPtr[k]];//Lik =/ Lkk

            //storing Lik = Lki
            for(int j4 = rowPtr[k] + 1; j4 < rowPtr[k+1]; ++j4){
                if(col[j4] == i){
                    lltvalues_[j4] = lltvalues_[j1];
                    break;
                }
                else if(col[j4] > i){
                    break;
                }
            }

            //Lik * Lik
            sumD -= lltvalues_[j1] * lltvalues_[j1];
        }

        if(sumD < 0){
            //to safeguard against negative sqr root
            //fallback to original diag value
            sumD = A.values()[rowPtr[i]];
        }
        lltvalues_[rowPtr[i]] = std::sqrt(sumD);
    }
}

void ic0Precond::apply
(
    std::vector<double>& z,
    const std::vector<double>& r
) const
{
    const auto& rowPtr = matPtr_->pattern()->rowPtr();
    const auto& col = matPtr_->pattern()->column();

    int n = r.size();
    
    //IC(0) Implementation
    //to solve: (LL^T)z = r
    //Let L^Tz = y -> solve Ly = r (forward elimination)
    //Then solve L^T*z = y (Backward substitution)
    
    for(int i = 0; i < n; ++i){//forward elimination
        double sum = r[i];
        for(int j = rowPtr[i] + 1; j < rowPtr[i+1]; ++j){

            if(col[j] > i) break;

            sum -= lltvalues_[j]*z[col[j]];
        }
        z[i] = sum/lltvalues_[rowPtr[i]];
    }

    for(int i = n-1; i >=0; --i){//backward substitution 
        double sum = z[i];
        for(int j = rowPtr[i] + 1; j < rowPtr[i+1]; ++j){

            if(col[j] < i) continue;

            sum -= lltvalues_[j] * z[col[j]];
        }
        z[i] = sum/lltvalues_[rowPtr[i]];
    }
}