#pragma once

#include "Applications/Applications.h"
#include "Grid/FvMesh.h"
#include "FiniteVolumeFields/ScalarField.h"

class Laplacian : public Applications
{
private:
    FvMesh mesh_;
    ScalarField phi_;
    double gamma_;

public:
    //Destructor
    virtual ~Laplacian() = default;

    //Constructor
    Laplacian();

    //initialize the scalar field
    void createField();

    //Prime Driver
    void run() override;
};