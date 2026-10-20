#include "Algebra/Solvers/PCG.h"
#include "Utils/VectorMath.h"

#include <cmath>

void PCG::solve
(
    const CsrMatrix& A,
    const std::vector<double>& b,
    std::vector<double>& x
) const
{
    //setting up M^-1 ~ A^-1
    precond_->setup(A);

    int n = b.size();
    
    std::vector<double> r(n, 0.0);//residual vector
    std::vector<double> z(n, 0.0);//preconditioned residual
    std::vector<double> p(n, 0.0);//search direction
    std::vector<double> Ax(n, 0.0);
    std::vector<double> Ap(n, 0.0);

    A.SpMV(x, Ax);
    VectorMath::sub(b, Ax, r);
    precond_->apply(z, r);
    p = z;

    //checking the initial residual
    double res = VectorMath::dot(r, r);
    res = std::sqrt(res);
    double res0 = res;
    if(res0 < tolerance_) return;

    double rz = VectorMath::dot(r, z);
    int iter = 0;

    while
    (
        iter < maxIter_ &&
        res > tolerance_ &&
        res/res0 > relTol_
    )
    {
        A.SpMV(p, Ap);

        double pAp = VectorMath::dot(p, Ap);
        double alpha = rz/pAp;

        //update x and residual
        for(int i = 0; i < n; ++i){
            x[i] = x[i] + alpha * p[i];
            r[i] = r[i] - alpha * Ap[i];
        }

        precond_->apply(z, r);

        double curr_rz = VectorMath::dot(r, z);
        
        double beta = curr_rz/rz;

        //update search direction
        for(int i = 0; i < n; ++i){
            p[i] = z[i] + beta*p[i];
        }

        rz = curr_rz;

        res = VectorMath::dot(r, r);
        res = std::sqrt(res);
        iter++;
    }
}