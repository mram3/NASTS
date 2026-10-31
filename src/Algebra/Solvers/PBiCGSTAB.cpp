#include "Algebra/Solvers/PBiCGSTAB.h"
#include "Utils/VectorMath.h"

#include <cmath>
#include <iostream>

SolverPerf PBiCGSTAB::solve
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
        zr_.assign(n, 0.0);//Bi residual vector
        p_.assign(n, 0.0);//search direction
        zp_.assign(n, 0.0);//preconditioned search direction
        s_.assign(n, 0.0);
        zs_.assign(n, 0.0);
        t_.assign(n, 0.0);
        Ax_.assign(n, 0.0);
        v_.assign(n, 0.0);//v = Ap~
    }
    
    //brief PBiCGStab (Barrett R. - Templates for solution of Linear Systems,
    //                 Building block for iterative methods)

    //First iteration is performed out of loop
    //r(0) = b - A*x(0)
    A.SpMV(x, Ax_);
    VectorMath::sub(b, Ax_, r_);
    //choosing r~ = r(0)
    zr_ = r_;

    //p1 = r0 (initial search direction taken to be r0)
    p_ = r_;

    //rho0 = r~^T * r(0) 
    double currRho = VectorMath::dot(zr_, r_);

    //inital residual
    double res = VectorMath::L2Norm(r_);
    double res0 = res;
    if(res0 < tolerance_) return SolverPerf{res0, res, 0};

    int iter = 0;
    while
    (
        iter < maxIter_ &&
        res > tolerance_ &&
        res > relTol_ * res0
    ){
        iter++;

        //preconditioning search direction
        precond_->apply(zp_, p_);

        //v = A*p~ (storing just for convinience)
        A.SpMV(zp_, v_);

        //step length alpha = r~^T * r / r~^T * A*p = rho / r~^T * v
        double zrv = VectorMath::dot(zr_, v_);
        if(std::abs(zrv) < 1e-35){
            std::cout << "PBiCGStab Method failed (alpha breakdown)\n";
            return SolverPerf{res0, res, iter};
        }
        double alpha = currRho/zrv;

        //s = r - alpha * v
        s_ = r_;
        VectorMath::daxpy(-alpha, v_, s_);

        //checking the norm of s
        double res_s = VectorMath::L2Norm(s_);
        if(res_s < relTol_*res0 || res_s < tolerance_){
            //if it's small enough, updating x and exiting
            VectorMath::daxpy(alpha, zp_, x);
            return SolverPerf{res0, res_s, iter};
        }

        //preconditioning s
        precond_->apply(zs_, s_);

        //t = A * preconditioned s
        A.SpMV(zs_, t_);

        //omega = t^T * s / t^T * t
        double tt = VectorMath::dot(t_, t_);

        if(tt < 1e-35){
            return SolverPerf{res0, res, iter};
        }
        double omega = VectorMath::dot(t_, s_)/tt;

        //for continuation it is necessary omega =/ 0
        if(std::isnan(omega) || std::abs(omega) < 1e-35){
            std::cout << "PBiCGStab method fails (omega breakdown)\n";
            return SolverPerf{res0, res, iter};
        }

        for(int i = 0; i < n; ++i){
            //updating x
            //x = x + alpha * p~ + omega * s~
            x[i] += alpha * zp_[i] + omega * zs_[i];
            
            //updating the residual 
            //r = s - omega * t
            r_[i] = s_[i] - omega * t_[i];
        }

        //calculate the Norm
        res = VectorMath::L2Norm(r_);
        //no need to update next direction if converged
        if(res < relTol_ * res0 || res < tolerance_){
            return SolverPerf{res0, res, iter};
        } 
        
        //to update the search direction for next iteration
        double prevRho = currRho;
        currRho = VectorMath::dot(zr_, r_);

        if(std::abs(currRho) < 1e-35){
            std::cout << "PBiCGStab Method failed (rho breakdown)\n";
            return SolverPerf{res0, res, iter};
        }

        //beta = rho(i)/rho(i-1) * (alpha/omega)
        double beta = (currRho/prevRho) * (alpha/omega);

        for(int i = 0; i < n; ++i){
            //update search direction
            //p(i) = r(i-1) + beta(i-1)[p(i-1) - omega(i-1)*v(i-1)]
            p_[i] = r_[i] + beta * (p_[i] - omega * v_[i]);
        }
    }

    return

    SolverPerf
    {
        res0, res, iter
    };
}