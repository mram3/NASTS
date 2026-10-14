#include "Grid/Grid.h"
#include "Matrix/CsrAddress.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <map>
#include <utility>
#include <cmath>

void Grid::readGmshAscii2
(
    std::ifstream& mshfile
)
{
    std::string line;

    std::vector<int> tempNodes;
    tempNodes.reserve(4);

    while(getline(mshfile, line)){

        if(line == "$MeshFormat"){
            getline(mshfile, line);
            std::stringstream ss(line);

            double version;
            ss >> version;

            if(version != 2.2) std::cout << "Unsupported .msh version";
        }

        else if(line == "$PhysicalNames"){
            getline(mshfile, line);
            int numNames = std::stoi(line);

            for(int i = 0; i < numNames; ++i){
                getline(mshfile, line);
                std::stringstream ss(line);

                int dim, tag; std::string name;

                ss >> dim >> tag >> name;

                //removing the " " Gmsh puts around physical names
                name.erase(std::remove(name.begin(), name.end(), '\"'), name.end());

                physicalName[tag] = name;
            }
        }

        else if(line == "$Nodes"){
            getline(mshfile, line);
            int numNodes = std::stoi(line);
            nodes.reserve(numNodes);

            for(int i = 0; i < numNodes; ++i){
                getline(mshfile, line);
                std::stringstream ss(line);

                int nodeId; double x, y, z;

                ss >> nodeId >> x >> y >> z;

                nodes.emplace_back(--nodeId, x, y, z);
            }
        }

        else if (line == "$Elements"){
            getline(mshfile, line);
            int numElements = std::stoi(line);
            elements.reserve(numElements);

            for(int i = 0; i < numElements; ++i){
                getline(mshfile, line);
                std::stringstream ss(line);

                int elementId, typeNum, tagCount, physicalTag, tag2;
                
                ss >> elementId >> typeNum >> tagCount >>physicalTag >> tag2;

                ElementType elementType = static_cast<ElementType>(typeNum);

                switch (elementType){
                    case ElementType::Line:{
                        int n1, n2;
                        ss >> n1 >> n2;

                        tempNodes = {--n1, --n2};
                        break;
                    }
                    case ElementType::Triangle:{
                        int n1, n2, n3;
                        ss >> n1 >> n2 >> n3;

                        tempNodes = {--n1, --n2, --n3};
                        break;
                    }
                    case ElementType::Quad:{
                        int n1, n2, n3, n4;
                        ss >> n1 >> n2 >> n3 >> n4;

                        tempNodes = {--n1, --n2, --n3, --n4};
                        break;
                    }
                }
                
                elements.emplace_back(elementId, elementType, physicalTag, tempNodes);
                tempNodes.clear();//is not required, just for my readability
            }
        }
    }
}

void Grid::buildMeshDetails
(

)
{
    faces.clear();
    cells.clear();

    //to map what pair of nodes constitute which face
    std::map<std::pair<int, int>, int> faceNodeMap;

    for(auto const &element : elements){

        if(element.elementType == ElementType::Line) continue;
        //just adding all the faces of all cells here, boundary faces are marked seperately..
        //.. below in an another loop

        Cell presentCell;

        presentCell.cellID = cells.size();
        
        int numVertices = element.nodeIDs.size();// 3 if triangle, 4 if quad etc.,

        presentCell.faceIDs.reserve(numVertices);

        for(int i = 0; i < numVertices; ++i){
            int n1 = element.nodeIDs[i];
            int n2 = element.nodeIDs[(i+1) % numVertices];//extracting two nodes at a time

            auto edge = std::minmax(n1, n2);//sorting the edge nodes

            //searching the map if the edge is there already
            if(faceNodeMap.find(edge) == faceNodeMap.end()){
                Face newface;//if not add a new face 

                int faceID = faces.size(); //temporary identifier for face

                newface.nodeA  = edge.first;
                newface.nodeB  = edge.second;

                newface.ownerCellID = presentCell.cellID;

                faces.push_back(newface);

                presentCell.faceIDs.push_back(faceID);

                faceNodeMap[edge] = faceID;
            }

            else{

                int oldfaceID = faceNodeMap[edge];//if edge is already present..

                faces[oldfaceID].neighborCellID = presentCell.cellID;//..update the rightside

                presentCell.faceIDs.push_back(oldfaceID);
            }
        }

        //storing the nodeIDs within Cell for easier geometry calculation
        presentCell.nodeIDs = element.nodeIDs;
        cells.push_back(presentCell);
    }

    for(const auto& element : elements){

        //marking the physical tag of boundary layer faces so it is easy to assign..
        //..boundary conditions during matrix assembly
        if(element.elementType == ElementType::Line){
            int n1 = element.nodeIDs[0];
            int n2 = element.nodeIDs[1];

            auto edge = std::minmax(n1, n2);

            auto id = faceNodeMap.find(edge);
            int faceID = id->second;

            faces[faceID].physicalTag = element.physicalTag;
        }
    }

    //sorting the faces array to have internal faces at front (prevents if clause during assembly)

    //partitioning internal faces to front
    auto split = std::partition(faces.begin(), faces.end(), [](const Face& face){
        return face.neighborCellID != -1;
    });

    numInternalFaces = std::distance(faces.begin(), split);

    //now since the faceIDs are displaced randomly to bring internal faces at front we
    //have to reassign those IDs to cell array
    std::vector<int> newID(faces.size());
    for(int i = 0; i < faces.size(); ++i){
        //now we retrive that temporary id in the map
        auto edge = std::minmax(faces[i].nodeA, faces[i].nodeB);
        int tempIdx = faceNodeMap[edge];
        newID[tempIdx] = i; //in place of that temp id we store the current id
    }

    for(auto& cell : cells){
        int i = 0;
        for(auto tempIdx : cell.faceIDs){
            cell.faceIDs[i++] = newID[tempIdx];
        }
    }
}

void Grid::evaluateGeometry
(

)
{
    //Cell Geometry Calculations
    for(auto& cell : cells){
        double xSum = 0, ySum = 0;
        int numNodes = cell.nodeIDs.size();

        //to calculate the center of each cell
        for(const auto& nodeID : cell.nodeIDs){
            xSum += nodes[nodeID].x;
            ySum += nodes[nodeID].y;
        }

        cell.cell_center_x = xSum / numNodes;
        cell.cell_center_y = ySum / numNodes;

        //volume
        double tempVol = 0;
        for(int i = 0; i < numNodes; ++i){
            int node1 = cell.nodeIDs[i];
            int node2 = cell.nodeIDs[(i+1) % numNodes];

            double x1 = nodes[node1].x;
            double x2 = nodes[node2].x;
            double y1 = nodes[node1].y;
            double y2 = nodes[node2].y;

            tempVol += x1*y2 - x2*y1;
        }

        cell.cell_volume = 0.5 * std::abs(tempVol);
    }

    //face Geometry calculations
    for(auto& face : faces){

        //to calculate center of a face
        double x1 = nodes[face.nodeA].x;
        double y1 = nodes[face.nodeA].y;
        double x2 = nodes[face.nodeB].x;
        double y2 = nodes[face.nodeB].y;

        face.face_center_x = 0.5 * (x1 + x2);
        face.face_center_y = 0.5 * (y1 + y2);

        //face area
        double dx = x1 - x2; 
        double dy = y1 - y2;
        face.face_area = std::sqrt(dx*dx + dy*dy);

        //outward facing normal
        double nx = - dy / face.face_area;
        double ny =   dx / face.face_area;

        double xLeft = cells[face.ownerCellID].cell_center_x;
        double yLeft = cells[face.ownerCellID].cell_center_y;

        face.xOf = face.face_center_x - xLeft;//vector from face center to cell center
        face.yOf = face.face_center_y - yLeft;

        double dot_prod = face.xOf*nx + face.yOf*ny;//dot product bw normal and celltoface vector

        if(dot_prod < 0) nx = -nx, ny = -ny;

        face.face_normal_x = nx;
        face.face_normal_y = ny;

        face.face_area_x = nx*face.face_area;
        face.face_area_y = ny*face.face_area;

        if(face.physicalTag == -1){
            face.xON = cells[face.neighborCellID].cell_center_x - xLeft;
            face.yON = cells[face.neighborCellID].cell_center_y - yLeft;
        }
    }
}

Grid::Grid
(
    std::ifstream& mshfile
)
{
    readGmshAscii2(mshfile);

    buildMeshDetails();

    evaluateGeometry();

    pattern_ = std::make_unique<CsrAddress>(*this);
}

Grid::~Grid() = default;