#pragma once

#include "FieldSpace.h"
#include "Type.h"

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

    //constructor
    ScalarField(std::string fieldName_, const Grid& mesh_): FieldSpace(fieldName_, mesh_) 
    {
        int numCells = mesh.getCells().size();

        scalarVariable.assign(numCells, 0.0);
        oldScalarVariable.assign(numCells, 0.0);
    };
    
    //helper function when i want to store old values of scalar quantities

};