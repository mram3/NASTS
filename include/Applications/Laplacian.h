#pragma once

#include "Applications/Applications.h"
#include "Grid/FvMesh.h"
#include "FiniteVolumeFields/ScalarField.h"

class Laplacian : public Applications
{
private:
    FvMesh fvmesh_;
    ScalarField phi_;
    double gamma_;

    /// @brief initialize the scalar field
    void createField() override;

    /// @brief check the user input file
    void check() const override;

public:
    //Destructor
    virtual ~Laplacian() = default;

    //Constructor
    Laplacian();

    //Prime Driver
    void run() override;

    const FvMesh& fvmesh() const {return fvmesh_; }
};