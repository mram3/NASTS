#pragma once

#include "Applications/Applications.h"
#include "Matrix/FvMatrix.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "Controls/Time.h"
#include "Controls/Output.h"

class UnsteadyLaplacian : public Applications
{
private:
    FvMesh fvmesh_;
    Time time_;
    Output output_;

    std::unique_ptr<ScalarField> phi_;
    double gamma_;
    double rho_;

    /// @brief initialize the scalar field
    void createField() override;

    /// @brief check the user input file
    void check() const override;

public:
    //Destructor
    virtual ~UnsteadyLaplacian() = default;

    //Constructor
    UnsteadyLaplacian();

    //The prime driver
    void run () override;
};