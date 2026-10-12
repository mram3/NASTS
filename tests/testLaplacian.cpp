/*
Compiling Instruction :g++ tests/testLaplacian.cpp src/*.cpp -Iinclude -std=c++17 -o test && ./test
*/
#include <iostream>
#include <cmath>

#include "Grid/Grid.h"
#include "Fields/ScalarField.h"
#include "Matrix/FvMatrix.h"
#include "Matrix/implicitAssembly.h"

using namespace std;

bool testLaplacian(const Grid& mesh){

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
                                                           2*x + 3*y};
        }
    }

    FvMatrix mat = implicitAssembly::laplacian(1.0, phi);

    mat.compressToCSR();

    std::vector<double> x = phi.scalarVariable;
    std::vector<double> lhs = mat.SpMV(x);

    std::vector<double> rhs = mat.getSource();

    double tolerance = 1e-4;

    for(auto const& cell : mesh.getCells()){
        if(std::abs(rhs[cell.cellID]-lhs[cell.cellID]) > tolerance){
            cout << "Failed at " << cell.cellID << endl;
            cout << lhs[cell.cellID] << "   |   " << rhs[cell.cellID] << endl;
            return false;
        }
    }

    return true;
}

int main(){

    //import the .msh file
    ifstream file("tests/assets/cavity.msh");

    if(!file.is_open()){
        cout << "\nCannot open the .msh file\n";
        return 1;
    }

    Grid mesh;

    mesh.readGmshAscii2(file);

    mesh.buildMeshDetails();

    mesh.evaluateGeometry();

    bool result  = testLaplacian(mesh);

    if(result){
        cout << "Laplacian Test Successfull\n";
    }
    else{
        cout << "Laplacian Test failed\n";
    }

    return 0;
}