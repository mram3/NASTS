/*
Compiling Instructions: g++ tests/testGrad.cpp src/*.cpp -Iinclude -std=c++17 -o test && ./test
*/
#include <iostream>
#include <fstream>
#include <cmath>

#include "Grid/Grid.h"
#include "Fields/ScalarField.h"
#include "Fields/VectorField.h"
#include "FiniteVolumeMethods/explicitCalc.h"

using namespace std;

bool testLinearGrad(const Grid& mesh){

    ScalarField phi("Phi", mesh);

    for(auto const& cell : mesh.getCells()){
        double x = cell.cell_center_x;
        double y = cell.cell_center_y;

        //test function: phi = 2*x + 3*y
        phi.scalarVariable[cell.cellID] = 2.0 * x + 3.0 * y;
    }

    //setting boundary conditions
    for(auto const& face : mesh.getFaces()){
        if(face.physicalTag != -1){
            double x = face.face_center_x;
            double y = face.face_center_y;

            phi.boundaryConditionMap[face.physicalTag] = {BoundaryConditionType::Dirichlet,
                                                           0};//some random value never accessed
        }
    }

    VectorField gradPhi = explicitCalc::grad(phi);

    double tolerance = 1e-8;

    for(const auto& cell : mesh.getCells()){
        double gradx = gradPhi.vectorVariable[cell.cellID].x;
        double grady = gradPhi.vectorVariable[cell.cellID].y;

        double dx = std::abs(gradx - 2.0);
        double dy = std::abs(grady - 3.0);
        if(dx > tolerance || dy > tolerance){
            cout << "Gradient is not matching at " << cell.cellID <<
                    ". Expected (2.0, 3.0) but got " << gradx <<", " << grady << endl;
            
            return false;
        }
    }

    return true;
}

int main(){

    //import the .msh file
    ifstream file("tests/assets/structured.msh");

    if(!file.is_open()){
        cout << "\nCannot open the .msh file\n";
        return 1;
    }

    Grid mesh(file);

    bool result = testLinearGrad(mesh);

    if(result){
        cout << "\nGradient Calculation Successfull\n";
    }
    else{
        cout << "\nGradient Calculation failed\n";
    }

    return 0;
}