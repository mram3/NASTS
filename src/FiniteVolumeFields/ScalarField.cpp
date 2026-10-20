#include "FiniteVolumeFields/ScalarField.h"

void ScalarField::create
(
    const nlohmann::json& fieldConfig
)
{
    int n = fvmesh_.cells().size();

    double iVal = fieldConfig.value("InitialCondition", 0.0);

    if(std::abs(iVal) > 1e-15){
        scalarVariable.assign(n, iVal);
        oldScalarVariable.assign(n, iVal);
    }

    auto bc = fieldConfig["BoundaryCondition"];

    for(const auto& item : bc.items()){
        std::string boundaryName = item.key();
        auto bcDetail = item.value();

        auto findTag = fvmesh_.physicalTag().find(boundaryName);

        if(findTag == fvmesh_.physicalTag().end()){
            throw std::runtime_error
            (
                "Can't find the boundary : " + boundaryName + " specified in input in .msh file"
            );
        }

        int tag = findTag->second;

        std::string bcType = bcDetail["Type"];
        double bcValue = bcDetail["Value"];

        if(bcType == "Dirichlet"){
            boundaryConditionMap[tag] = 
            {
                BoundaryConditionType::Dirichlet, 
                bcValue
            };
        }
        else if(bcType == "Neumann"){
            boundaryConditionMap[tag] = 
            {
                BoundaryConditionType::Neumann, 
                bcValue
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