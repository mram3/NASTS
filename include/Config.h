#pragma once

#include "nlohmann/json.hpp"

#include <fstream>
#include <string>
#include <stdexcept>

class Config
{
private:
    //the only container which stores all information provided in input json file
    nlohmann::json data_;

    //private constructor (Singleton)
    Config() {}

    //To check if the absolutely necessary information are in usr input before starting simulation
    void checkInput() const;
    
public:

    //Access to the Config
    static Config& access()
    {
        static Config inst;
        return inst;
    }

    //To read the user json file at the start of the program via access()
    void read(const std::string& filePath)
    {
        std::ifstream file(filePath);

        if(!file.is_open()){
            throw std::runtime_error("Input file not found: " + filePath);
        }
        data_ = nlohmann::json::parse(file);

        //Immediately validate the input file before moving on to heavy functions
        checkInput();
    }

    //getter to send user provided information to multiple classes via access()
    const nlohmann::json& data() const {return data_; }
};