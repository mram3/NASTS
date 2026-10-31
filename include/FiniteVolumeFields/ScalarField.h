#pragma once

#include "FieldSpace.h"
#include "Type.h"
#include "Config.h"
#include "nlohmann/json.hpp"

#include <vector>
#include <map>
#include <string>

class ScalarField : public FieldSpace
{
private:
    int size_;
    bool initTooOld_ = false;

public:

    //storage for physical variables like pressure, Temperature etc.,
    std::vector<double> scalarVariable;
    std::vector<double> oldScalarVariable;
    std::vector<double> tooOldScalarVariable;

    //mapping user provided bc tag and the bc type and value
    std::map<int, BoundaryCondition<double>> boundaryConditionMap;

    //construct from FvMesh and Field Name
    ScalarField
    (
        std::string fieldName, 
        const FvMesh& mesh
    )
    :   FieldSpace(fieldName, mesh),
        size_(mesh.cells().size())
    {
        scalarVariable.assign(size_, 0.0);

        std::string scheme = Config::access().data()["Time"]["Scheme"];

        if(scheme == "SOUE") initTooOld_ = true;
    };

    //creates scalar field from config
    void create
    (
        const nlohmann::json& fieldConfig
    );
    
    //helper function when i want to store old values of scalar quantities
    void storeField();
};