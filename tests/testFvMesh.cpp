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
        const auto& boundaries = fvmesh.fvboundaries();
        const auto& domains = fvmesh.fvdomains();

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
            boundaries.at("Lid") == 5 &&
            boundaries.at("Walls") == 6 &&
            domains.at("InternalField") == 7,
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