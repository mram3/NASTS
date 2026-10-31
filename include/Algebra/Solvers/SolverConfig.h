#pragma once

#include "Algebra/Solvers/LinearSolvers.h"
#include "Config.h"

#include <unordered_map>
#include <string>
#include <memory>

class SolverConfig
{
private:

    //map that stores the requested solver for a specific field
    mutable std::unordered_map
    <
        std::string,
        std::unique_ptr<LinearSolvers>

    > solvers_;

    //private constructor - SingleTon
    SolverConfig() {}

public:

    static SolverConfig& access(){
        static SolverConfig inst;
        return inst;
    }

    LinearSolvers& solver
    (
        std::string& fieldName
    )const 
    {
        auto settings = Config::access().data()["LinearSolvers"][fieldName];

        if(solvers_.find(fieldName) == solvers_.end()){
            solvers_[fieldName] = LinearSolvers::make(settings);

            return *solvers_[fieldName];
        }

        return *solvers_[fieldName];
    }
};