#include "Grid/FvMesh.h"
#include "Grid/MeshReader.h"
#include "Matrix/FvAddressing.h"

#include <map>
#include <algorithm>
#include <cmath>

FvMesh::FvMesh
(
    
)
: FvMesh(MeshReader::readGmshAscii22())
{
    //constructor 
}

FvMesh::FvMesh
(
    const Grid& mesh
)
{
    //constructor
    
    //FVM faces_ and cells are built from nodes and elements
    buildFvMesh(mesh);
    evalGeometry(mesh);

    //map phyicalName to tag
    for(const auto& pair : mesh.boundaries()){
        fvboundaries_[pair.second] = pair.first;
    }
    //map phyicalName to tag
    for(const auto& pair : mesh.domians()){
        fvdomains_[pair.second] = pair.first;
    }
}

void FvMesh::buildFvMesh
(
    const Grid& mesh
)
{
    faces_.clear();
    cells_.clear();

    const auto& elements = mesh.elements();

    //to map what pair of nodes constitute which face
    std::map<std::pair<int, int>, int> faceNodeMap;

    for(auto const &element : elements){

        if(element.elementType == ElementType::Line) continue;
        //just adding all the faces of all cells here, boundary faces are marked seperately..
        //.. below in an another loop

        Cell presentCell;

        presentCell.cellID = cells_.size();
        
        int numVertices = element.nodeIDs.size();// 3 if triangle, 4 if quad etc.,

        presentCell.faceIDs.reserve(numVertices);

        for(int i = 0; i < numVertices; ++i){
            int n1 = element.nodeIDs[i];
            int n2 = element.nodeIDs[(i+1) % numVertices];//extracting two nodes at a time

            auto edge = std::minmax(n1, n2);//sorting the edge nodes

            //searching the map if the edge is there already
            if(faceNodeMap.find(edge) == faceNodeMap.end()){
                Face newface;//if not add a new face 

                int faceID = faces_.size(); //temporary identifier for face

                newface.nodeA  = edge.first;
                newface.nodeB  = edge.second;

                newface.ownerCellID = presentCell.cellID;

                faces_.push_back(newface);

                presentCell.faceIDs.push_back(faceID);

                faceNodeMap[edge] = faceID;
            }

            else{

                int oldfaceID = faceNodeMap[edge];//if edge is already present..

                faces_[oldfaceID].neighborCellID = presentCell.cellID;//..update the rightside

                presentCell.faceIDs.push_back(oldfaceID);
            }
        }

        //storing the nodeIDs within Cell for easier geometry calculation
        presentCell.nodeIDs = element.nodeIDs;
        cells_.push_back(presentCell);
    }

    for(const auto& element : elements){

        //marking the physical tag of boundary faces so it is easy to assign..
        //..boundary conditions during matrix assembly
        if(element.elementType == ElementType::Line){
            int n1 = element.nodeIDs[0];
            int n2 = element.nodeIDs[1];

            auto edge = std::minmax(n1, n2);

            auto id = faceNodeMap.find(edge);
            int faceID = id->second;

            faces_[faceID].physicalTag = element.physicalTag;
        }
    }

    //sorting the faces array to have internal faces at front 
    //(prevents if clause during assembly)

    //partitioning internal faces to front
    auto split = std::partition(faces_.begin(), faces_.end(), [](const Face& face){
        return face.neighborCellID != -1;
    });

    numInternalFaces_ = std::distance(faces_.begin(), split);

    //now since the faceIDs are displaced randomly to bring internal faces at front we
    //have to reassign those IDs to cell array
    std::vector<int> newID(faces_.size());
    for(int i = 0; i < faces_.size(); ++i){
        //now we retrive that temporary id in the map
        auto edge = std::minmax(faces_[i].nodeA, faces_[i].nodeB);
        int tempIdx = faceNodeMap[edge];
        newID[tempIdx] = i; //in place of that temp id we store the current id
    }

    for(auto& cell : cells_){
        int i = 0;
        for(auto tempIdx : cell.faceIDs){
            cell.faceIDs[i++] = newID[tempIdx];
        }
    }
}

void FvMesh::evalGeometry
(
    const Grid& mesh
)
{
    nodes_ = mesh.nodes();
    const auto& elements = mesh.elements();
    //Cell Geometry Calculations
    for(auto& cell : cells_){
        double xSum = 0, ySum = 0;
        int numNodes = cell.nodeIDs.size();

        //to calculate the center of each cell
        for(const auto& nodeID : cell.nodeIDs){
            xSum += nodes_[nodeID].x;
            ySum += nodes_[nodeID].y;
        }

        cell.cell_center_x = xSum / numNodes;
        cell.cell_center_y = ySum / numNodes;

        //volume
        double tempVol = 0;
        for(int i = 0; i < numNodes; ++i){
            int node1 = cell.nodeIDs[i];
            int node2 = cell.nodeIDs[(i+1) % numNodes];

            double x1 = nodes_[node1].x;
            double x2 = nodes_[node2].x;
            double y1 = nodes_[node1].y;
            double y2 = nodes_[node2].y;

            tempVol += x1*y2 - x2*y1;
        }

        cell.cell_volume = 0.5 * std::abs(tempVol);
    }

    //face Geometry calculations
    for(auto& face : faces_){

        //to calculate center of a face
        double x1 = nodes_[face.nodeA].x;
        double y1 = nodes_[face.nodeA].y;
        double x2 = nodes_[face.nodeB].x;
        double y2 = nodes_[face.nodeB].y;

        face.face_center_x = 0.5 * (x1 + x2);
        face.face_center_y = 0.5 * (y1 + y2);

        //face area
        double dx = x1 - x2; 
        double dy = y1 - y2;
        face.face_area = std::sqrt(dx*dx + dy*dy);

        //outward facing normal
        double nx = - dy / face.face_area;
        double ny =   dx / face.face_area;

        double xLeft = cells_[face.ownerCellID].cell_center_x;
        double yLeft = cells_[face.ownerCellID].cell_center_y;

        face.xOf = face.face_center_x - xLeft;//vector from face center to cell center
        face.yOf = face.face_center_y - yLeft;

        double dot_prod = face.xOf*nx + face.yOf*ny;//dot product bw normal & celltoface vector

        if(dot_prod < 0) nx = -nx, ny = -ny;

        face.face_normal_x = nx;
        face.face_normal_y = ny;

        face.face_area_x = nx*face.face_area;
        face.face_area_y = ny*face.face_area;

        if(face.physicalTag == -1){
            face.xON = cells_[face.neighborCellID].cell_center_x - xLeft;
            face.yON = cells_[face.neighborCellID].cell_center_y - yLeft;
        }
    }
}

const std::shared_ptr<SparseAddress> FvMesh::fvAddress
(

) const
{
    if(!fvAddress_){
        fvAddress_ = std::make_shared<SparseAddress> 
                    (
                        FvAddressing::buildFvAddress(*this)
                    );
    }

    return fvAddress_;
}

const std::shared_ptr<SparseAddress> FvMesh::fvDiagAddress
(

) const
{
    if(!fvDiagAddress_){
        fvDiagAddress_ = std::make_shared<SparseAddress> 
                        (
                            FvAddressing::buildFvDiagAddress(*this)
                        );
    }

    return fvDiagAddress_;
}