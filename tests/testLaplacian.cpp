#include "Tests.h"
#include "Grid/FvMesh.h"
#include "Grid/MeshReader.h"
#include "Matrix/FvMatrix.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "FiniteVolumeMethods/implicitAssembly/implicitAssembly.h"
#include "Algebra/Solvers/LinearSolvers.h"
#include "nlohmann/json.hpp"

#include <iostream>
#include <string>
#include <cmath>

int main(){
    std::cout << "--Laplacian Test--\n";

    try
    {
        //solve a field using Dirichlet BC
        Grid mesh = MeshReader::readGmshAscii22("tests/assets/square.msh");
        FvMesh fvmesh(mesh);

        ScalarField T1("T1", fvmesh);

        std::string field1 = R"(
        {
            "InitialCondition" : 0.0,
            "BoundaryCondition" : {
                "LeftAndRight" : {"Type" : "Dirichlet", "Value" : 500},
                "TopAndBottom"  : {"Type" : "Dirichlet", "Value" : 300.0} 
            }
        })";

        nlohmann::json fieldConfig1 = nlohmann::json::parse(field1);
        T1.create(fieldConfig1);

        FvMatrix fvMat1 = implicitAssembly::laplacian(1.0, T1);

        fvMat1 == 0.0;

        std::string solver1 = R"(
        {
            "Solver" : "PCG",
            "Preconditioner" : "Jacobi",
            "MaxIter" : 1000,
            "Tolerance" : 1e-6,
            "RelTol" : 0.0
        })";

        nlohmann::json solverConfig1 = nlohmann::json::parse(solver1);
        auto linearSolver1 = LinearSolvers::make(solverConfig1);

        linearSolver1->solve(fvMat1, fvMat1.source(), T1.scalarVariable);

        //test dirichlet BC
        const auto& result1 = T1.scalarVariable;

        EXPECT_TRUE
        (
            std::abs(result1[12] - 400.0) < 1e-10,
            "Laplacian Assembly and Solve Method with Dirichlet BC"
        );

        //solve a field with Neumann BC
        ScalarField T2("T2", fvmesh);
        std::string field2 = R"(
        {
            "InitialCondition" : 0.0,
            "BoundaryCondition" : {
                "LeftAndRight" : {"Type" : "Dirichlet", "Value" : 500},
                "TopAndBottom"  : {"Type" : "Neumann", "Value" : 0.0} 
            }
        })";

        nlohmann::json fieldConfig2 = nlohmann::json::parse(field2);
        T2.create(fieldConfig2);

        FvMatrix fvMat2 = implicitAssembly::laplacian(1.0, T2);

        fvMat2 == 0.0;

        std::string solver2 = R"(
        {
            "Solver" : "PCG",
            "Preconditioner" : "Jacobi",
            "MaxIter" : 1000,
            "Tolerance" : 1e-6,
            "RelTol" : 0.0
        })";

        nlohmann::json solverConfig2 = nlohmann::json::parse(solver2);
        auto linearSolver2 = LinearSolvers::make(solverConfig2);

        linearSolver2->solve(fvMat2, fvMat2.source(), T2.scalarVariable);

        //test dirichlet BC
        const auto& result2 = T2.scalarVariable;
        
        EXPECT_TRUE
        (
            std::abs(result2[12] - 500.0) < 1e-10,
            "Laplacian Assembly and Solve Method with Neumann BC"
        );
    }
    catch(const std::exception& e)
    {
        std::cout << "\033[31m[ CRASH ]\033[0m " << e.what() << "\n";
        return 1;
    }
    
    return Tests::summary();
}