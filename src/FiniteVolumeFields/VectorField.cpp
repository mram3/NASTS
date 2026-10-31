#include "FiniteVolumeFields/VectorField.h"
#include "Type.h"

#include <string>

void VectorField::createVelocityField
(
    const nlohmann::json& fieldConfig
)
{
    int n = fvmesh_.cells().size();

    auto init = fieldConfig["InitialCondition"];

    double iValx = init.value("u", 0.0);
    double iValy = init.value("v", 0.0);
    TwoVector iVal = {iValx, iValy};

    vectorVariable.assign(n, iVal);
    oldVectorVariable.assign(n, iVal);

    auto bc = fieldConfig["BoundaryCondition"];

    for(const auto& item : bc.items()){
        std::string boundaryName = item.key();
        auto bcDetail = item.value();

        int tag = fvmesh_.fvboundaries().at(boundaryName);

        std::string bcType = bcDetail["Type"];
        auto bcValue = bcDetail["Value"];

        TwoVector boundary = {bcValue["u"], bcValue["v"]};

        if(bcType == "Dirichlet"){
            boundaryConditionMap[tag] = 
            {
                BoundaryConditionType::Dirichlet, 
                boundary
            };
        }

        else if(bcType == "Neumann"){
            boundaryConditionMap[tag] = 
            {
                BoundaryConditionType::Neumann, 
                boundary
            };
        }

        else{
            throw std::runtime_error
            (
                "Unknown boundary condition specified : " + bcType
            );
        }
    }
}