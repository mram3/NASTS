#include "FiniteVolumeMethods/implicitAssembly/ConvectionSchemes/ConvectionSchemes.h"
#include "FiniteVolumeMethods/implicitAssembly/ConvectionSchemes/Upwind.h"
#include "FiniteVolumeMethods/implicitAssembly/ConvectionSchemes/DeferredCorrection.h"
#include "FiniteVolumeMethods/surfaceInterpolation/vanLeer.h"
#include "FiniteVolumeMethods/surfaceInterpolation/MUSCL.h"
#include "FiniteVolumeMethods/implicitAssembly/implicitAssembly.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "FiniteVolumeFields/VectorField.h"
#include "Config.h"

#include <iostream>
#include <memory>

namespace implicitAssembly
{

FvMatrix div
(
    SurfaceScalarField& flux,
    ScalarField& phi
)
{
    auto convConfig = Config::access().data()["Convection"]["Scheme"];

    static auto divPtr = ConvectionScheme::make
    (
        convConfig
    );

    return
    divPtr->div
    (
        flux,
        phi
    );
}

std::unique_ptr<ConvectionScheme> ConvectionScheme::make
(
    const std::string& scheme
)
{
    if(scheme == "Upwind"){

        return
        std::make_unique<Upwind>();
    }
    else if(scheme == "vanLeer"){

        auto interpScheme = std::make_unique<fvinterp::vanLeer>();

        return 
        std::make_unique<DeferredCorrection>
        (
            std::move(interpScheme)
        );
    }
    else if(scheme == "MUSCL"){
        auto interpScheme = std::make_unique<fvinterp::MUSCL>();

        return 
        std::make_unique<DeferredCorrection>
        (
            std::move(interpScheme)
        );
    }
    else{
        throw std::runtime_error
        (
            "Unknown Convection Scheme: " + scheme
        );
    }
}
}