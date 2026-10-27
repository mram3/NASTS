#pragma once

#include "Matrix/FvMatrix.h"
#include "FiniteVolumeFields/ScalarField.h"
#include <memory>
#include <string>

namespace implicitAssembly
{

/// @brief Factory class to choose the time stepping scheme once at start up
class DdtSchemes
{
public:
    //Destructor
    virtual ~DdtSchemes() = default;

    /// @brief Prime Matrix assembly method defined as pure virtual function here
    /// @param rho - Density
    /// @param phi - Scalar Quantity being solved
    /// @return - Assembled FvMatrix (Diagonal matrix)
    virtual FvMatrix Ddt
    (
        double rho,
        ScalarField& phi
    ) 
    = 0;

    /// @brief Factory method to initiate the requested scheme class
    /// @param schemeName - Name of the scheme
    /// @return - Pointer to requested scheme class
    static std::unique_ptr<DdtSchemes> make
    (
        const std::string& schemeName
    );
};

}