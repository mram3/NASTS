#pragma once

#include "Matrix/CsrMatrix.h"
#include "Algebra/Preconditioners/Preconditioners.h"
#include "nlohmann/json.hpp"
#include "Type.h"

#include <memory>
#include <vector>

class LinearSolvers
{
protected:
    std::unique_ptr<Preconditioners> precond_;
    double maxIter_;
    double relTol_;
    double tolerance_;

public:

    //Destructor
    virtual ~LinearSolvers() = default;

    //Constructor
    LinearSolvers
    (
        std::unique_ptr<Preconditioners> precond,
        double maxIter,
        double relTol,
        double tolerance
    )
    :   precond_(std::move(precond)), //move the ownership to Preconditioners class
        maxIter_(maxIter),
        relTol_(relTol),
        tolerance_(tolerance)
    {}

    //The prime solver method which solves our A*x = b iteratively
    virtual SolverPerf solve
    (
        const CsrMatrix& A,
        const std::vector<double>& b,
        std::vector<double>& x
    ) = 0;

    //Solvers Factory
    static std::unique_ptr<LinearSolvers> make
    (
        const nlohmann::json& file
    );
};