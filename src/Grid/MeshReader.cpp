#include "Grid/MeshReader.h"
#include "Config.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace MeshReader
{

Grid readGmshAscii22
(
    const std::string& filePath
)
{
    std::cout << "Reading the mesh file " + filePath + "\n";

    std::ifstream mshFile(filePath);

    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::map<int, std::string> boundaries;
    std::map<int, std::string> domains;

    std::string line;

    std::vector<int> tempNodes;
    tempNodes.reserve(4);//assuming quad cells

    while(getline(mshFile, line)){

        if(line == "$MeshFormat"){
            getline(mshFile, line);
            std::stringstream ss(line);

            double version;
            ss >> version;

            if(version != 2.2){
                throw std::runtime_error("Error: unsupported mesh file type: " 
                                          + std::to_string(version));
            }
        }

        else if(line == "$PhysicalNames"){
            getline(mshFile, line);
            int numNames = std::stoi(line);

            for(int i = 0; i < numNames; ++i){
                getline(mshFile, line);
                std::stringstream ss(line);

                int dim, tag; std::string name;

                ss >> dim >> tag;

                //To handle spaces between names 
                std::getline(ss >> std::ws, name);

                //removing the " " Gmsh puts around physical names
                name.erase(std::remove(name.begin(), name.end(), '\"'), name.end());

                if(dim == 1) boundaries[tag] = name;

                else if(dim ==2) domains[tag] = name;
            }
        }

        else if(line == "$Nodes"){
            getline(mshFile, line);
            int numNodes = std::stoi(line);
            nodes.reserve(numNodes);

            for(int i = 0; i < numNodes; ++i){
                getline(mshFile, line);
                std::stringstream ss(line);

                int nodeId; double x, y, z;

                ss >> nodeId >> x >> y >> z;

                nodes.emplace_back(--nodeId, x, y, z);
            }
        }

        else if (line == "$Elements"){
            getline(mshFile, line);
            int numElements = std::stoi(line);
            elements.reserve(numElements);

            for(int i = 0; i < numElements; ++i){
                getline(mshFile, line);
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

    mshFile.close();

    std::cout << "Mesh File read\n\n";

    //construct grid and return
    return
    Grid
    (
        std::move(nodes),
        std::move(elements),
        std::move(boundaries),
        std::move(domains)
    );
}

Grid readGmshAscii22
(
    
)
{
    std::string meshFile = Config::access().data()["Mesh"];

    return readGmshAscii22(meshFile);
}
}