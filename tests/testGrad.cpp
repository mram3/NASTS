#include "Tests.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "FiniteVolumeFields/VectorField.h"
#include "Grid/FvMesh.h"
#include "Grid/Grid.h"
#include "Grid/MeshReader.h"
#include "Type.h"
#include "FiniteVolumeMethods/explicitCalc/explicitCalc.h"

#include <iostream>
#include <string>

bool testResult(VectorField &result);

int main(){
    std::cout << "--Explicit Gradient Calculation Test--\n";

    try
    {
        Grid mesh = MeshReader::readGmshAscii22("tests/assets/fourBoundaries.msh");

        FvMesh fvmesh(mesh);

        ScalarField phi("T", fvmesh);

        std::string field = R"(
        {
            "InitialCondition" : 500.0,
            "BoundaryCondition" : {
                "Left" : {"Type" : "Dirichlet", "Value" : 0.0},
                "Right"  : {"Type" : "Dirichlet", "Value" : 200.0},
                "Top" : {"Type" : "Neumann", "Value" : 0.0},
                "Bottom" : {"Type" : "Neumann", "Value" : 0.0}  
            }
        })";

        nlohmann::json fieldConfig = nlohmann::json::parse(field);
        phi.create(fieldConfig);

        const auto& cells = fvmesh.cells();
        for(const auto& cell : cells){
            phi.scalarVariable[cell.cellID] = 200*cell.cell_center_x;
        }

        VectorField result = explicitCalc::grad(phi);

        EXPECT_TRUE
        (
            testResult(result),
            "Explicit Gradient Calculation"
        );
    }
    catch(const std::exception& e)
    {
        std::cout << "\033[31m[ CRASH ]\033[0m " << e.what() << "\n";
        return 1;
    }
    
    return Tests::summary();
}

bool testResult
(
    VectorField& result
)
{
    for(int i = 0; i < result.vectorVariable.size(); ++i){
        if(
            std::abs(result.vectorVariable[i].x - 200.0) > 1e-6 ||
            std::abs(result.vectorVariable[i].y - 0.0) > 1e-6
        ){
            std::cout << "Fails at " << i << " Calculated Value : " 
                      << result.vectorVariable[i].x << " " << result.vectorVariable[i].y;
            return false;
        }
    }

    return true;
}