#include "Config.h"

#include <fstream>
#include <string>
#include <stdexcept>

void Config::checkInput
(

) const
{
    if(!data_.contains("Application"))
    {
        throw std::runtime_error
        (
            "CONFIG. ERROR : \"Application\" name is not specified\n"
        );
    }

    if(!data_["Application"].is_string())
    {
        throw std::runtime_error
        (
            "CONFIG ERROR: Application name has to be wrapped in \"quotes\"\n"
        );
    }

    if(!data_.contains("Mesh"))
    {
        throw std::runtime_error
        (
            "CONFIG. ERROR : \"Mesh\" Path is not specified\n"
        );
    }

    if(!data_["Mesh"].is_string())
    {
        throw std::runtime_error
        (
            "CONFIG ERROR: Mesh Path has to be wrapped in \"quotes\"\n"
        );
    }

    std::string meshPath = data_["Mesh"];
    std::ifstream meshFile(meshPath);
    if(!meshFile.is_open())
    {
        throw std::runtime_error
        (
            "CONFIG. ERROR : Cannot open the \"Mesh\" file at " + meshPath + "\n"
        );
    }
    meshFile.close();

    if(!data_.contains("Fields"))
    {
        throw std::runtime_error
        (
            "CONFIG. ERROR : \"Fields\" section is not specified\n"
        );
    }
    
    if(!data_["Fields"].is_object())
    {
        throw std::runtime_error
        (
            "CONFIG ERROR: The \"Fields\" section should be a dictionary wrapped in "
            "braces {.....}\n"
        );
    }

    if(data_["Fields"].empty())
    {
        throw std::runtime_error
        (
            R"("CONFIG. ERROR : "Fields" section is empty. 
                Specify Initial and Boundary Conditions.
                )"
        );
    }

    if(!data_.contains("LinearSolvers"))
    {
        throw std::runtime_error
        (
            "CONFIG. ERROR : \"LinearSolvers\" section is not specified\n"
        );
    }

    if(!data_["LinearSolvers"].is_object())
    {
        throw std::runtime_error
        (
            "CONFIG ERROR: The \"LinearSolvers\" section should be a dictionary wrapped in "
            "braces {.....}\n"
        );
    }

    for(auto& [name, settings] : data_["LinearSolvers"].items())
    {
        if(!data_["Fields"].contains(name))
        {
            throw std::runtime_error
            (
                "CONFIG ERROR: Missing \"" + name + "\" defined in \"LinearSolvers\"\n in "
                "\"Fields\" Block"
            );
        }

        if(!settings.is_object())
        {
            throw std::runtime_error
            (
                "CONFIG ERROR: The " + name + " field section inside "
                "\"LinearSolvers\" section should be a dictionary wrapped in "
                "braces {.....}\n"
            );
        }

        if(!settings.contains("Solver"))
        {
            throw std::runtime_error
            (
                "CONFIG ERROR: No \"Solver\" is specified for the field : " + name + "\n"
            );
        }
    }
}