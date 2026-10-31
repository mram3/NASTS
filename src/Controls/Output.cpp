#include "Controls/Output.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <algorithm>

void Output::writeVtk
(
    bool force
) const
{
    //write only at the requested write interval or if forced to write
    if(!force && !time_.write()){ return;}

    const auto& cells = fvmesh_.cells();
    const auto& nodes = fvmesh_.nodes();

    std::string outDir = "VTK";
    if(!std::filesystem::exists(outDir)){
        std::filesystem::create_directory(outDir);
    }

    int frameIdx = time_.isSteadySolve() ? time_.currIter() 
                                         : time_.currStep();

    std::stringstream fileName;
    fileName <<  outDir << "/output_" << std::setfill('0') 
             << std::setw(4) << frameIdx << ".vtk";

    std::ofstream vtkFile(fileName.str());
    if(!vtkFile.is_open()){
        std::cerr << "Unable to open VTK file " << fileName.str() << " for writing\n";
        return;
    }

    //VTK HEADER
    vtkFile << "# vtk DataFile Version 3.0\n";
    if(time_.isSteadySolve()){
        vtkFile << "NASTS CFD Output at iteration " << time_.currIter() << "\n";
    }
    else{
        vtkFile << "NASTS CFD Output at time " << time_.currTime() << "\n";
    }
    vtkFile << "ASCII\n";
    vtkFile << "DATASET UNSTRUCTURED_GRID\n";

    //POINTS (NODES)
    vtkFile << "POINTS " << nodes.size() << " double\n";
    for(const auto& node : nodes){
        vtkFile << node.x << " " << node.y << " " << node.z << "\n";
    }

    //CELLS (CONNECTIVITY)
    int totCellDetails = 0;
    for(const auto& cell : cells){
        // 1 for type + N nodes per cell
        totCellDetails += 1 + cell.nodeIDs.size();
    }

    vtkFile << "\nCELLS " << cells.size() << " " << totCellDetails << "\n";
    for(const auto& cell : cells){
        //no of nodes in this cell
        vtkFile << cell.nodeIDs.size();
        for(int nodeID : cell.nodeIDs){
            vtkFile << " " << nodeID;
        }
        vtkFile << "\n";
    }

    //CELL TYPES
    vtkFile << "\nCELL_TYPES " << cells.size() << "\n";
    for(const auto& cell : cells){
        int numNodes = cell.nodeIDs.size();
        int cellType = 0;

        if(numNodes == 3){
            cellType = 5; //Triangle in VTK
        }
        else if(numNodes == 4){
            cellType = 9; //QUAD in VTK
        }
        else{
            cellType = 7; //Generic Polygon in VTK
        }
        
        vtkFile << cellType << "\n";
    }

    //CELL DATA
    vtkFile << "\nCELL_DATA " << cells.size() << "\n";
    
    //writing all scalars
    for(const auto* phi : scalarFields_){
        vtkFile << "SCALARS " << phi->fieldName() << " double 1\n";
        vtkFile << "LOOKUP_TABLE default\n";

        for(const auto& cell : cells){
            vtkFile << phi->scalarVariable[cell.cellID] << "\n";
        }
    }

    //write all vectors
    for(const auto* U : vectorFields_){
        vtkFile << "VECTORS " << U->fieldName() << " double\n";

        for(const auto& cell : cells){
            vtkFile << U->vectorVariable[cell.cellID].x << " " 
                    << U->vectorVariable[cell.cellID].y << " " << "0.0\n";
        }
    }

    vtkFile.close();
}