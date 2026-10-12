#pragma once

#include "Grid/Grid.h"
#include "Type.h"

#include <string>

class FieldSpace
{
protected:
    std::string fieldName;
    const Grid& mesh;

public:

    //Destructor
    virtual ~FieldSpace() = default;

    //Constrctor
    FieldSpace(std::string fieldName_, const Grid& mesh_): fieldName(fieldName_), mesh(mesh_) {}

    //getter functions 
    std::string getFieldName() const { return fieldName; }
    const Grid& getMesh() const { return mesh; }
};