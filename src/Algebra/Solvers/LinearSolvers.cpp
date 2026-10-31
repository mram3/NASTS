#include "Algebra/Solvers/LinearSolvers.h"
#include "Algebra/Solvers/PCG.h"
#include "Algebra/Solvers/PBiCGSTAB.h"
#include "Algebra/Preconditioners/GaussFamPrecond.h"
#include "Algebra/Preconditioners/IncompleteFamPrecond.h"
#include "Algebra/Preconditioners/noPrecond.h"
#include "Config.h"

#include <memory>
#include <string>
#include <stdexcept>
#include <iostream>

std::unique_ptr<LinearSolvers> LinearSolvers::make
(
    const nlohmann::json& file
)
{
    //Solver settings
    int maxIter = file.value("MaxIter", 1000);
    double relTol  = file.value("RelTol", 0.0);
    double tolerance = file.value("Tolerance", 1e-6);

    //preconditioner type
    std::unique_ptr<Preconditioners> precond;
    //check for user specified preconditioner. By default 'No Preconditioner' is applied
    std::string pcType = file.value("Preconditioner", "NONE");

    if(pcType == "Jacobi"){
        precond = std::make_unique<jacobiPrecond>();
    }
    else if(pcType == "SGS"){
        precond  = std::make_unique<sgsPrecond>();
    }
    else if(pcType == "ILU0"){
        precond  = std::make_unique<ilu0Precond>();
    }
    else if(pcType == "IC0"){
        precond  = std::make_unique<ic0Precond>();
    }
    else{
        precond = std::make_unique<noPrecond>();
    }

    //solver type
    std::string solver = file["Solver"];

    if(solver == "PCG"){
        return 
        
        std::make_unique<PCG>
        (
            std::move(precond),
            maxIter,
            relTol,
            tolerance
        );
    }

    else if(solver == "PBiCGStab"){
        return

        std::make_unique<PBiCGSTAB>
        (
            std::move(precond),
            maxIter,
            relTol,
            tolerance
        );
    }

    throw std::runtime_error
    (
        "Unknown solver type specified : " + solver 
    );
}