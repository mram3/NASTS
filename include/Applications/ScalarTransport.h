#pragma once

#include "Applications/Applications.h"
#include "Matrix/FvMatrix.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "FiniteVolumeFields/SurfaceScalarField.h"
#include "FiniteVolumeFields/VectorField.h"
#include "Controls/Time.h"
#include "Controls/Output.h"

#include <memory>

class ScalarTransport : public Applications
{
private:
    FvMesh fvmesh_;
    Time time_;
    Output output_;

    std::unique_ptr<ScalarField> phi_;
    std::unique_ptr<SurfaceScalarField> flux_;
    std::unique_ptr<VectorField> U_;
    
    double gamma_;
    double rho_;

    /// @brief initialize the scalar field
    void createField() override;

    /// @brief check the user input file
    void check() const override;

public:
    //Destructor
    virtual ~ScalarTransport() = default;

    //Constructor
    ScalarTransport();

    //The prime driver
    void run () override;
};