#include "Algebra/Solvers/PCG.h"
#include "Utils/VectorMath.h"

#include <cmath>

SolverPerf PCG::solve
(
    const CsrMatrix& A,
    const std::vector<double>& b,
    std::vector<double>& x
)
{
    //setting up M^-1 ~ A^-1
    precond_->setup(A);

    int n = b.size();
    
    if(r_.empty()){
        r_.assign(n, 0.0);//residual vector
        z_.assign(n, 0.0);//preconditioned residual
        p_.assign(n, 0.0);//search direction
        Ax_.assign(n, 0.0);
        Ap_.assign(n, 0.0);
    }

    A.SpMV(x, Ax_);
    VectorMath::sub(b, Ax_, r_);
    precond_->apply(z_, r_);
    p_ = z_;

    //checking the initial residual
    double res = VectorMath::L2Norm(r_);
    double res0 = res;
    if(res0 < tolerance_) return SolverPerf{res0, res, 0};

    double rz = VectorMath::dot(r_, z_);
    int iter = 0;

    while
    (
        iter < maxIter_ &&
        res > tolerance_ &&
        res > relTol_ * res0
    )
    {
        A.SpMV(p_, Ap_);

        double pAp = VectorMath::dot(p_, Ap_);
        double alpha = rz/pAp;

        //update x and residual
        for(int i = 0; i < n; ++i){
            x[i] = x[i] + alpha * p_[i];
            r_[i] = r_[i] - alpha * Ap_[i];
        }

        precond_->apply(z_, r_);

        double curr_rz = VectorMath::dot(r_, z_);
        
        double beta = curr_rz/rz;

        //update search direction
        for(int i = 0; i < n; ++i){
            p_[i] = z_[i] + beta*p_[i];
        }

        rz = curr_rz;

        res = VectorMath::dot(r_, r_);
        res = std::sqrt(res);
        iter++;
    }

    return 
    SolverPerf
    {
        res0, res, iter
    };
}