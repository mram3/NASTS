#pragma once

#include <memory>
#include <string>

class Applications
{
public:
    //Destructor
    virtual ~Applications() = default;

    virtual void check() const = 0;

    virtual void createField() = 0;

    virtual void run() = 0;

    static std::unique_ptr<Applications> make
    (
        const std::string& appType
    );
};