#pragma once

#include "Grid/FvMesh.h"
#include "Type.h"

#include <string>

class FieldSpace
{
protected:
    std::string fieldName_;
    const FvMesh& fvmesh_; //referance to the mesh of the field

public:

    //Destructor
    virtual ~FieldSpace() = default;

    //Constrctor
    FieldSpace
    (
        std::string fieldName, 
        const FvMesh& mesh
    )
    :   fieldName_(fieldName),
        fvmesh_(mesh) 
    {}

    //getter functions 
    std::string fieldName() const { return fieldName_; }
    const FvMesh& fvmesh() const { return fvmesh_; }
};