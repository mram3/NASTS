#pragma once

#include "FieldSpace.h"
#include "Type.h"

#include <map>

class VectorField : public FieldSpace
{
public:
    //storage containers to store vector variables of 2 components
    std::vector<TwoVector> vectorVariable;
    std::vector<TwoVector> oldVectorVariable;

    //mapping the user defined bc tag to bc type and value
    std::map<int, BoundaryCondition<TwoVector>> boundaryConditionMap;

    //Constructor
    VectorField(std::string fieldName_, const Grid& mesh) : FieldSpace(fieldName_, mesh)
    {
        int numCells = mesh.getCells().size();

        vectorVariable.assign(numCells, {0.0, 0.0});
        oldVectorVariable.assign(numCells, {0.0, 0.0});
    }
    
};