#include "Applications/Applications.h"
#include "Applications/Laplacian.h"
#include "Applications/UnsteadyLaplacian.h"

#include <memory>
#include <stdexcept>

std::unique_ptr<Applications> Applications::make
(
    const std::string& appType
)
{
    if(appType == "Laplacian"){
        return

        std::make_unique<Laplacian>();
    }

    else if(appType == "UnsteadyLaplacian"){
        return

        std::make_unique<UnsteadyLaplacian>();
    }

    throw std::runtime_error("Unknown Application requested : " + appType);
}