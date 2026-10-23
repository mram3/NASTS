#include "Applications/Laplacian.h"
#include "Applications/UnsteadyLaplacian.h"
#include "Grid/FvMesh.h"
#include "Config.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
void checkGamma(const nlohmann::json& data)
{
    if(!data.contains("Gamma"))
    {
        throw std::runtime_error
            (
                "APP. ERROR: Missing \"Gamma\" value\n"
            );
    }

    if(!data["Gamma"].is_number())
    {
        throw std::runtime_error
        (
            "App. ERROR: Gamma value should NOT be wrapped in \"quotes\"\n"
        );
    }
}

void checkRho(const nlohmann::json& data)
{
    if(!data.contains("Rho"))
    {
        throw std::runtime_error
            (
                "APP. ERROR: Missing \"Rho\" value\n"
            );
    }

    if(!data["Rho"].is_number())
    {
        throw std::runtime_error
        (
            "App. ERROR: Rho value should NOT be wrapped in \"quotes\"\n"
        );
    }
}

void checkSource(const nlohmann::json& data)
{
    if(!data.contains("Source"))
    {
        std::cout << "App. WARNING: \"Source\" term is not specified. Default Value is ZERO\n";
    }
    else
    {
        if(!data["Source"].is_number())
        {
            throw std::runtime_error
            (
                "App. ERROR: \"Source\" value should NOT be wrapped in \"quotes\"\n"
            );
        }
    }
}

void checkNOC(const nlohmann::json& data)
{
    if(!data.contains("nNonOrthogonalCorrections"))
    {
        std::cout << "App. WARNING: \"nNonOrthogonalCorrections\" term is not specified. "  
                     "Default Value is ZERO\n";
    }
    else
    {
        if(!data["nNonOrthogonalCorrections"].is_number())
        {
            throw std::runtime_error
            (
                "App. ERROR: \"nNonOrthogonalCorrections\" value should NOT be"
                "wrapped in \"quotes\"\n"
            );
        }
    }
}

void checkBC(const nlohmann::json& data, const FvMesh& fvmesh)
{
    std::string msh = data["Mesh"];

    for(const auto& [name, settings] : data["Fields"].items())
    {
        if(!settings.contains("InitialCondition"))
        {
            throw std::runtime_error
            (
                "APP. ERROR: Missing \"InitialCondition\" block for " + name + 
                " field\n"
            );
        }

        if(!settings.contains("BoundaryCondition"))
        {
            throw std::runtime_error
            (
                "APP. ERROR: Missing \"BoundaryCondition\" block for " + name + 
                " field\n"
            );
        }

        for(const auto& [mshBName, tag] : fvmesh.fvboundaries())
        {
            if(!settings["BoundaryCondition"].contains(mshBName))
            {
                throw std::runtime_error
                (
                    "APP. ERROR: \"" + mshBName + "\" boundary mentioned in " + msh + 
                    " file is missing in input file\n"
                );
            }
        }

        for(const auto& boundary : settings["BoundaryCondition"].items())
        {

            std::string bName = boundary.key();
            auto bSettings = boundary.value();

            if(fvmesh.fvboundaries().find(bName) == fvmesh.fvboundaries().end())
            {
                throw std::runtime_error
                (
                    "APP. ERROR: \"" + bName + "\" boundary mentioned in input file " 
                    "is missing in " + msh + " file\n"
                );
            }

            if(!bSettings.contains("Type"))
            {
                throw std::runtime_error
                (
                "APP. ERROR: Missing \"Type\" block for " + bName + 
                " boundary\n"
                "Supported Types are \"Dirichlet\" and \"Neumann\"\n"
                );
            }

            if(!bSettings.contains("Value"))
            {
                throw std::runtime_error
                (
                "APP. ERROR: Missing \"Value\" block for " + bName + 
                " boundary\n"
                );
            }

            std::string bcType = bSettings["Type"];
            if(bcType != "Dirichlet" && bcType != "Neumann")
            {
                throw std::runtime_error
                (
                "APP. ERROR: UnKnown \"Type\" block for " + bName + 
                " boundary\n"
                "Supported Types are \"Dirichlet\" and \"Neumann\"\n"
                );
            }
        }
    }
}

void checkTime(const nlohmann::json& data)
{
    if(!data.contains("Time"))
    {
        throw std::runtime_error
        (
            "APP. ERROR: Missing \"Time\" block\n"
        );
    }

    if(!data["Time"].is_object())
    {
        throw std::runtime_error
        (
            "CONFIG ERROR: The \"Time\" section should be a dictionary wrapped in "
            "braces {.....}\n"
        );
    }

    if(!data["Time"].contains("Scheme"))
    {
        throw std::runtime_error
        (
            "APP. ERROR: Missing \"Scheme\" inside \"Time\" block\n"
        );
    }

    if(!data["Time"]["Scheme"].is_string())
    {
        throw std::runtime_error
        (
            "CONFIG ERROR: The \"Scheme\" inside \"Time\" block should be a wrapped in "
            "\"quotes\"\n"
        );
    }

    if(!data["Time"].contains("StartTime"))
    {
        throw std::runtime_error
        (
            "APP. ERROR: Missing \"StartTime\" inside \"Time\" block\n"
        );
    }

    if(!data["Time"]["StartTime"].is_number())
    {
        throw std::runtime_error
        (
            "CONFIG ERROR: The \"StartTime\" inside \"Time\" block should NOT be a wrapped in "
            "\"quotes\"\n"
        );
    }

    if(!data["Time"].contains("EndTime"))
    {
        throw std::runtime_error
        (
            "APP. ERROR: Missing \"EndTime\" inside \"Time\" block\n"
        );
    }

    if(!data["Time"]["EndTime"].is_number())
    {
        throw std::runtime_error
        (
            "CONFIG ERROR: The \"EndTime\" inside \"Time\" block should NOT be a wrapped in "
            "\"quotes\"\n"
        );
    }

    if(!data["Time"].contains("dt"))
    {
        throw std::runtime_error
        (
            "APP. ERROR: Missing \"dt\" inside \"Time\" block\n"
        );
    }

    if(!data["Time"]["dt"].is_number())
    {
        throw std::runtime_error
        (
            "CONFIG ERROR: The \"dt\" inside \"Time\" block should NOT be a wrapped in "
            "\"quotes\"\n"
        );
    }
}

void checkSteadyState(const nlohmann::json& data, bool ss)
{
    if(ss)
    {
    if(!data["Time"]["SteadyState"].contains("Tolerance"))
    {
        std::cout << "APP. WARNING: Steady State tolerance not specified. "
                     "Default value is 1e-6\n";
    }
    }
}
}

void Laplacian::check
(

) const
{
    const auto& data = Config::access().data();

    checkGamma(data);
    checkSource(data);
    checkNOC(data);
    checkBC(data, fvmesh_);
}

void UnsteadyLaplacian::check
(

) const
{
    const auto& data = Config::access().data();

    checkGamma(data);
    checkSource(data);
    checkNOC(data);
    checkTime(data);
    checkBC(data, fvmesh_);
}