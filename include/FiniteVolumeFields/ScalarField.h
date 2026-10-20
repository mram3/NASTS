#pragma once

#include "FieldSpace.h"
#include "Type.h"
#include "nlohmann/json.hpp"

#include <vector>
#include <map>

class ScalarField : public FieldSpace
{
public:

    //storage for physical variables like pressure, Temperature etc.,
    std::vector<double> scalarVariable;
    std::vector<double> oldScalarVariable;

    //mapping user provided bc tag and the bc type and value
    std::map<int, BoundaryCondition<double>> boundaryConditionMap;

    //construct from FvMesh and Field Name
    ScalarField
    (
        std::string fieldName, 
        const FvMesh& mesh
    )
        : FieldSpace(fieldName, mesh) 
    {
        int numCells = mesh.cells().size();

        scalarVariable.assign(numCells, 0.0);
        oldScalarVariable.assign(numCells, 0.0);
    };

    //creates scalar field from config
    void create
    (
        const nlohmann::json& fieldConfig
    );
    
    //helper function when i want to store old values of scalar quantities

};