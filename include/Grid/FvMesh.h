#pragma once

#include "Grid/Grid.h"

#include <vector>
#include <memory>
#include <map>

struct Face
{
    int nodeA, nodeB; //what defines the face
    int ownerCellID = -1;  
    int neighborCellID = -1; //-1 would mean boundary face
    int physicalTag = -1; //-1 would mean internal face
    
    double face_center_x, face_center_y;

    double face_area; // *length in case of 2D (just written area for cfd scheme readability)
    double face_area_x, face_area_y; //components of face area

    double face_normal_x, face_normal_y;
    double xOf, yOf; //distance between owner cell center to face center

    double xON = 0.0, yON = 0.0;//distance between owner and neighbor(defalut to 0 for boundaries) 
};

struct Cell
{
    int cellID;
    std::vector<int> faceIDs; //faceIDs surrounding the cell
    std::vector<int> nodeIDs; //nodeIDs on the corners of cell

    double cell_center_x, cell_center_y;
    double cell_volume; // *area in case of 2D 
};

class SparseAddress;//link between the mesh and matrix

class FvMesh
{
private:

    std::vector<Face> faces_;
    std::vector<Cell> cells_;
    size_t numInternalFaces_ = 0;

    void buildFvMesh(const Grid& mesh);
    void evalGeometry(const Grid& mesh);

    mutable std::shared_ptr<SparseAddress> fvAddress_;

    mutable std::shared_ptr<SparseAddress> fvDiagAddress_;

    //to initialize boundary condition during field setup
    std::map<std::string, int> fvboundaries_;
    std::map<std::string, int> fvdomains_;

    //to write vtk files during post processing
    std::vector<Node> nodes_;

public:

    //construct from configuration
    FvMesh();

    //construct from mesh (for test purposes)
    explicit FvMesh(const Grid& mesh);

    //getter for field assembly/initialization
    const std::map<std::string, int>& fvboundaries() const {return fvboundaries_; }
    const std::map<std::string, int>& fvdomains() const {return fvdomains_; }

    //getters for FvMatrix assembly
    const std::vector<Face>& faces() const { return faces_; }
    const std::vector<Cell>& cells() const { return cells_; }
    const int numInternalfaces() const { return numInternalFaces_; }

    //Delayed Initialization (created only when a computational method function asks for it)
    const std::shared_ptr<SparseAddress> fvAddress() const;

    const std::shared_ptr<SparseAddress> fvDiagAddress() const;
};