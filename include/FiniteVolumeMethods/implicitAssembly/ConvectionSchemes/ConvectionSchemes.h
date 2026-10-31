#pragma once

#include "Matrix/FvMatrix.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "FiniteVolumeFields/VectorField.h"
#include "FiniteVolumeFields/SurfaceScalarField.h"

#include <memory>
#include <string>

namespace implicitAssembly
{

class ConvectionScheme
{
public:
    //Destructor
    virtual ~ConvectionScheme() = default;

    /// @brief Prime Matrix assembly method declared as pure virtual function here
    /// @param rho - Density of the fluid
    /// @param phi - Scalar field
    /// @param U - Velocity field
    /// @return - Assembled FV Matrix 
    virtual FvMatrix div
    (
        SurfaceScalarField& flux,
        ScalarField& phi
    ) const
    = 0;

    static std::unique_ptr<ConvectionScheme> make
    (
        const std::string& convectionConfig
    );
};

}