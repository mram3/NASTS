#include "FiniteVolumeMethods/implicitAssembly/implicitAssembly.h"
#include "FiniteVolumeMethods/implicitAssembly/DdtSchemes/DdtSchemes.h"
#include "FiniteVolumeMethods/implicitAssembly/DdtSchemes/EulerDdt.h"

#include "Config.h"

#include <memory>

namespace implicitAssembly

{

FvMatrix Ddt
(
    double rho,
    ScalarField& phi
)
{
    std::string scheme = Config::access().data()["Time"]["Scheme"];

    auto ddtSchemePtr = DdtSchemes::make(scheme);

    return ddtSchemePtr->Ddt(rho, phi);
}

std::unique_ptr<DdtSchemes> DdtSchemes::make
(
    const std::string& schemeName
)
{
    if(schemeName == "ImplicitEuler"){
        return

        std::make_unique<EulerDdt>();
    }

    else{
        throw 

        std::runtime_error("Unknown Time Stepping Scheme\n" + schemeName +
                           "Available schemes are \"ImplicitEuler\"");
    }
}
}