#include "Tests.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "Grid/FvMesh.h"
#include "Grid/Grid.h"
#include "Grid/MeshReader.h"
#include "Type.h"

#include <iostream>
#include <string>
#include <cmath>

int main(){
    
    std::cout << "--Scalar Fields Test--\n";

    try{

        //setting up the field
        Grid mesh = MeshReader::readGmshAscii22("tests/assets/cavity.msh");

        FvMesh fvmesh(mesh);

        ScalarField phi("T", fvmesh);

        std::string field = R"(
        {
            "InitialCondition" : 500.0,
            "BoundaryCondition" : {
                "Lid" : {"Type" : "Dirichlet", "Value" : 200},
                "Walls"  : {"Type" : "Neumann", "Value" : 0.0} 
            }
        })";

        nlohmann::json fieldConfig = nlohmann::json::parse(field);

        phi.create(fieldConfig);

        //tests
        EXPECT_TRUE
        (
            std::abs(phi.scalarVariable[0] - 500.0) < 1e-15,
            "Successful Scalar Field Initialization"
        );

        EXPECT_TRUE
        (
            std::abs(phi.oldScalarVariable[0] - 500.0) < 1e-15,
            "Successful Old Scalar Field Initialization"
        );

        EXPECT_TRUE
        (
            std::abs(phi.boundaryConditionMap[5].value - 200) < 1e-15 &&
            std::abs(phi.boundaryConditionMap[6].value - 0.0) < 1e-15,
            "Boundary Condition Value"
        );

        EXPECT_TRUE
        (
            phi.boundaryConditionMap[5].bc == BoundaryConditionType::Dirichlet &&
            phi.boundaryConditionMap[6].bc == BoundaryConditionType::Neumann,
            "Boundary Condition Type"
        );

    }catch(const std::exception& e){
        std::cout << "\033[31m[ CRASH ]\033[0m " << e.what() << "\n";
        return 1;
    }

    return Tests::summary();
}