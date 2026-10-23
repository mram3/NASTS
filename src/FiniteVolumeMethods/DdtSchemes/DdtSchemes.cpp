#include "FiniteVolumeMethods/implicitAssembly.h"
#include "Config.h"

namespace implicitAssembly
{
FvMatrix Ddt
(
    double rho,
    ScalarField& U
)
{
    std::string scheme = Config::access().data()["Time"]["Scheme"];

    if(scheme == "ImplicitEuler") return EulerDdt(rho, U);

    else{
        throw std::runtime_error (
            "Unknown Time Stepping Scheme : " + scheme
        );
    }
}
}