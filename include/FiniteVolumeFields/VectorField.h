#pragma once

#include "FieldSpace.h"
#include "Type.h"

#include <vector>
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
    VectorField
    (
        std::string fieldName, 
        const FvMesh& mesh
    ) 
        : FieldSpace(fieldName, mesh)
    {
        int numCells = mesh.cells().size();

        vectorVariable.assign(numCells, {0.0, 0.0});
        oldVectorVariable.assign(numCells, {0.0, 0.0});
    }
    
};