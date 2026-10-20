#pragma once

#include <memory>
#include <string>

class Applications
{
public:
    //Destructor
    virtual ~Applications() = default;

    virtual void run() = 0;

    static std::unique_ptr<Applications> make
    (
        const std::string& appType
    );
};