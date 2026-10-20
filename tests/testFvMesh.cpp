#include "Tests.h"
#include "Grid/Grid.h"
#include "Grid/MeshReader.h"
#include "Grid/FvMesh.h"

#include <iostream>

bool internalFaceTest(const std::vector<Face>& faces, int n);

double vol(const std::vector<Cell>& cells);

bool faceArea(const std::vector<Cell>& cells, const std::vector<Face>& faces);

bool faceReassignment(const std::vector<Cell>& cells, const std::vector<Face>& faces);

int main(){
    std::cout << "--FvMesh Test--\n";

    try{

        Grid mesh = MeshReader::readGmshAscii22("tests/assets/cavity.msh");

        FvMesh fvmesh(mesh);

        const auto& faces = fvmesh.faces();
        const auto& cells = fvmesh.cells();
        const auto& numInFace = fvmesh.numInternalfaces();
        const auto& physicalTag = fvmesh.physicalTag();

        EXPECT_EQ
        (
            faces.size(),
            109,
            "No. of Faces"
        );

        EXPECT_EQ
        (
            cells.size(),
            66,
            "No. of Cells"
        );

        EXPECT_EQ
        (
            numInFace,
            89,
            "No. of Internal Faces"
        );

        EXPECT_TRUE
        (
            internalFaceTest(faces, numInFace),
            "Internal Faces are arranged at front, followed by Boundary faces"
        );

        EXPECT_TRUE
        (
            faceReassignment(cells, faces),
            "Face Reassignment"
        );

        EXPECT_EQ
        (
            vol(cells),
            1.0,
            "Total Volume"
        );

        EXPECT_TRUE
        (
            faceArea(cells, faces),
            "Area and Normal Calculations"
        );

        EXPECT_TRUE
        (
            physicalTag.at("Lid") == 5 &&
            physicalTag.at("Walls") == 6 &&
            physicalTag.at("InternalField") == 7,
            "Physical Tag"
        );

    } catch(const std::exception& e){
        std::cout << "\033[31m[ CRASH ]\033[0m " << e.what() << "\n";
        return 1;
    }

    return Tests::summary();
    
}

bool internalFaceTest(const std::vector<Face>& faces, int n){
    for(int i = 0; i < n; ++i){
        if(faces[i].neighborCellID == -1) return false;
    }

    for(int j = n; j < faces.size(); ++j){
        if(faces[j].neighborCellID != -1) return false;
    }

    return true;
}

bool faceReassignment
(
    const std::vector<Cell>& cells,
    const std::vector<Face>& faces
)
{
    for(const auto& cell : cells){
        //we go to each cell
        for(auto faceID : cell.faceIDs){
            //for each face that cell either has to be owner or neighbor
            if(
                faces[faceID].ownerCellID != cell.cellID &&
                faces[faceID].neighborCellID != cell.cellID
            ) return false;
        }
    }

    return true;
}

double vol(const std::vector<Cell>& cells){
    double sum = 0.0;
    for(const auto& cell: cells){
        sum += cell.cell_volume;
    }

    return sum;
}

bool faceArea
(
    const std::vector<Cell>& cells,
    const std::vector<Face>& faces
)
{
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

    if(std::abs(areaXSum) > 1e-15 || std::abs(areaYSum) > 1e-15) return false;

    return true;
}
/*
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
         << " is " << pressure << endl;*/

