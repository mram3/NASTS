#include "FiniteVolumeMethods/implicitAssembly/implicitAssembly.h"
#include "FiniteVolumeMethods/implicitAssembly/DdtSchemes/DdtSchemes.h"
#include "FiniteVolumeMethods/implicitAssembly/DdtSchemes/EulerDdt.h"
#include "FiniteVolumeMethods/implicitAssembly/DdtSchemes/SOUEDdt.h"
#include "FiniteVolumeMethods/implicitAssembly/DdtSchemes/SteadyStateDdt.h"

#include "Config.h"

#include <memory>
#include <iostream>

namespace implicitAssembly

{

FvMatrix Ddt
(
    double rho,
    ScalarField& phi
)
{
    static auto ddtSchemePtr = DdtSchemes::make
    (
        Config::access().data()["Time"]["Scheme"]
    );

    return ddtSchemePtr->Ddt(rho, phi);
}

std::unique_ptr<DdtSchemes> DdtSchemes::make
(
    const std::string& schemeName
)
{

    if(schemeName == "SteadyState"){
        return

        std::make_unique<SteadyStateDdt>();
    }
    else if(schemeName == "ImplicitEuler"){
        return

        std::make_unique<EulerDdt>();
    }
    else if(schemeName == "SOUE"){
        return

        std::make_unique<SOUEDdt>();
    }

    else{
        throw 

        std::runtime_error("Unknown Time Stepping Scheme\n" + schemeName +
                           "\nAvailable schemes are \"ImplicitEuler\", \"SOUE\" and "
                           "\"SteadyState\"\n");
    }
}
}