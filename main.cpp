/*
Compiling Instruction: g++ main.cpp src/*.cpp -Iinclude -std=c++17 -o main && ./main
*/

#include <iostream>
#include <fstream>

#include "Grid/Grid.h"
#include "Fields/FieldSpace.h"
#include "Fields/ScalarField.h"
#include "Fields/VectorField.h"

using namespace std;

int main(){

    ifstream file("cavity.msh");

    if(!file.is_open()){
        cout << "Cannot open the .msh file\n";
        return 1;
    }

    Grid mesh;

    mesh.readGmshAscii2(file);
    file.close();

    auto nodes = mesh.getNodes();
    auto elements = mesh.getElements();
    auto phyNames = mesh.getPhysicalName();

    //checking what are read
    cout << "--Mesh Summary--\n";

    cout << "Number of nodes read: " << nodes.size() << endl;
    cout << "Number of elements read: " << elements.size() << endl;
    cout << "Number of Physical names read: " << phyNames.size() << endl;

    //checking what are the Physical names read
    cout << "\n--Physical Names--\n";

    for(const auto& [tag, name] : phyNames){
        cout << tag << " " <<phyNames[tag] << endl;
    }

    //checking the first node
    cout << "\n--First Node--\n";
    
    if(nodes.size() != 0){
        cout << nodes[0].nodeID << " " << nodes[0].x
                                << " " << nodes[0].y
                                << " " << nodes[0].z << endl;
    }

    mesh.buildMeshDetails();

    auto& faces = mesh.getFaces();
    auto& cells = mesh.getCells();

    //checking if topology function counted cells and faces correctly
    cout << "\n--Topology Details--\n";
    cout << "Number of faces: " << faces.size() << endl;
    cout << "Number of cells: " << cells.size() << endl;

    int boundary = 0;
    int invalidBoundary = 0;
    int invalidInteriorFace = 0;
    for(auto& face : faces){
        if(face.physicalTag != -1){
            ++boundary;

            if(face.neighborCellID != -1) ++invalidBoundary;
        }
        else{
            if(face.neighborCellID == -1) ++invalidInteriorFace;
        }
    }

    //boundary + interior cells should equal total elements in .msh file
    elements.size() == boundary + cells.size() ? cout << "\nElements are matching\n" 
                                               : cout <<"\nElements are not matching\n";

    cout << "Invalid Boundary Count : " << invalidBoundary << endl;
    cout << "Invalid Interior Faces Count : " << invalidInteriorFace << endl;

    mesh.evaluateGeometry();

    cout << "\n--Geometry Details--\n";

    //checking the sum of volumes are equal to total volume 
    double vol = 0;
    for(auto& cell : cells){
        vol  += cell.cell_volume;
    }

    cout << "Total volume: " << vol << endl;

    //the sum of dot products of outward normal and area of all faces of a polygon should be 0
    double areaXSum = 0, areaYSum = 0;
    for(auto& cell: cells){
        for(auto& faceID: cell.faceIDs){

            double nx = faces[faceID].face_normal_x;
            double ny = faces[faceID].face_normal_y;
            if(faces[faceID].neighborCellID == cell.cellID) nx = -nx, ny = -ny;
            areaXSum +=  nx * faces[faceID].face_area;
            areaYSum +=  ny * faces[faceID].face_area;
        }
    }

    cout << "Sum of x face areas: " << areaXSum << endl;
    cout << "Sum of y face areas: " << areaYSum << endl;

    //checking the Field class

    ScalarField p("Pressure", mesh);
    VectorField U("Velocity", mesh);

    for(const auto& face : faces){
        if(face.physicalTag == -1){
            p.boundaryConditionMap[-1] = {BoundaryConditionType::None, 0.0};
            U.boundaryConditionMap[-1] = {BoundaryConditionType::None, {0.5, 0.5}};
        }
        else if(face.physicalTag == 5){
            p.boundaryConditionMap[5] = {BoundaryConditionType::Neumann, 0.0};
            U.boundaryConditionMap[5] = {BoundaryConditionType::Dirichlet, {1.0, 0.0}};
        }
        else if(face.physicalTag == 6){
            p.boundaryConditionMap[6] = {BoundaryConditionType::Dirichlet, 0.0};
            U.boundaryConditionMap[6] = {BoundaryConditionType::Dirichlet, {0.0, 0.0}};
        }
    }

    cout << "\n--Field Details--\n";

    cout << "\nSuccessfully filled " << p.scalarVariable.size() << " cells with " 
         << p.getFieldName() << " field\n";

    cout << "\nSuccessfully filled " << U.vectorVariable.size() << " cells with " 
         << U.getFieldName() << " field\n";

    auto id = U.boundaryConditionMap.find(5);
    double u = id->second.value.x;
    double v = id->second.value.y;

    cout << "\nBoundary Condition for " << U.getFieldName() << " at " << phyNames[5]
         << " is " << u << " in x dirn and " << v << " in y dirn\n";

    auto it = p.boundaryConditionMap.find(5);
    double pressure = it->second.value;

    cout << "\nBoundary Condition for " << p.getFieldName() << " at " << phyNames[5]
         << " is " << pressure << endl;

    return 0;
}