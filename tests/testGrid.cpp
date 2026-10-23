/*
Compiling Instruction: g++ main.cpp src/*.cpp -Iinclude -std=c++17 -o main && ./main
*/
#include "Tests.h"
#include "Grid/Grid.h"
#include "Grid/MeshReader.h"

#include <iostream>
#include <fstream>

using namespace std;

int main(){

    std::cout << "--Grid Test--\n";

    try{
        Grid mesh = MeshReader::readGmshAscii22("tests/assets/cavity.msh");

        const auto& nodes = mesh.nodes();
        const auto& elements = mesh.elements();
        const auto boundaries = mesh.boundaries();
        const auto domains = mesh.domians();

        EXPECT_EQ
        (
            nodes.size(),
            44,
            "No. of Nodes"
        );

        EXPECT_EQ
        (
            elements.size(),
            86,
            "No. of Elements"
        );

        EXPECT_EQ
        (
            boundaries.size(),
            2,
            "No. of Physical Names"
        );

        EXPECT_EQ
        (
            domains.size(),
            1,
            "No. of Physical Names"
        );

        //checking if the physical names are mapped correctly
        EXPECT_TRUE
        (
            boundaries.at(5) == "Lid" &&
            boundaries.at(6) == "Walls" &&
            domains.at(7) == "InternalField",

            "Physical Name"
        );

        //checking the first node
        EXPECT_TRUE
        (
            nodes[0].nodeID == 0 &&
            std::abs(nodes[0].x) < 1e-15 &&
            std::abs(nodes[0].y) < 1e-15 &&
            std::abs(nodes[0].z) < 1e-15,
            
            "First Node"
        );
    }catch(const std::exception& e){
        std::cout << "\033[31m[ CRASH ]\033[0m " << e.what() << "\n";
        return 1;
    }

    return Tests::summary();
}