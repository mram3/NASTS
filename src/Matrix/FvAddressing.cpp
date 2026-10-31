#include "Matrix/FvAddressing.h"

#include <vector>
#include <algorithm>
#include <iostream>

namespace FvAddressing
{

SparseAddress buildFvAddress
(
    const FvMesh& mesh
)
{
    std::cout << "Building the fvAddress\n";

    const auto& cells = mesh.cells();
    const auto& faces = mesh.faces();

    size_t nCells = cells.size();
    size_t nInFace = mesh.numInternalfaces();

    //rowPtr is nothing but a count of non-zero elements before that row
    std::vector<int> rowPtr_(nCells + 1, 0);

    //so we lookup for number of non-zero contribution by each cell we have, by looking at how
    //many internal faces are shared by each cell and ofcourse a cell will have it's own
    //non-zero coefficient aP...

    std::vector<int> nnz(nCells, 1); //...so we start with 1

    for(size_t i = 0; i < nInFace; ++i){//boundary faces do not contribute to matrix A
        auto& face = faces[i];

        nnz[face.ownerCellID]++; 
        nnz[face.neighborCellID]++; 
    }

    //now that we have for each cell how many nnz it has, now we just progressively add it
    //to get nnz before each row

    rowPtr_[0] = 0; //has to be 
    for(size_t i = 0; i < nCells; ++i){
        rowPtr_[i+1] = rowPtr_[i] + nnz[i];
    }

    std::vector<int> column_(rowPtr_[nCells], 0);

    std::vector<int> tracker = rowPtr_; //to track within a row
    for(size_t i = 0; i < nCells; ++i){//just the cell's own column (i.e diag entry)
        column_[tracker[i]] = cells[i].cellID;
        tracker[i]++;
    }

    for(size_t i = 0; i < nInFace; ++i){
        const auto& face = faces[i];
        int oID = face.ownerCellID;
        int nID = face.neighborCellID;

        column_[tracker[oID]] = nID;//each face's owner cell row will have contribution from
                                    //neighbor cell, thus will have a column entry there
        column_[tracker[nID]] = oID;//and vice versa

        tracker[oID]++;
        tracker[nID]++;
    }
    
    //sorting the column's now keeping the diagonal entries first
    for(size_t i = 0; i < nCells; ++i){
        std::sort
        (
            column_.begin() + rowPtr_[i]+1, 
            column_.begin() + rowPtr_[i+1]
        );
    }

    std::cout << "FvAddress Built\n\n";

    return 
    SparseAddress
    (
        std::move(rowPtr_),
        std::move(column_)
    );
}

SparseAddress buildFvDiagAddress
(
    const FvMesh& fvmesh
)
{
    const auto& cells = fvmesh.cells();
    size_t n = cells.size();

    std::vector<int> rowPtr_(n+1, 0);
    std::vector<int> column_(n, 0);

    for(size_t i = 0; i < n; ++i){
        rowPtr_[i] = i;
        column_[i] = i;
    }
    rowPtr_[n] = n;

    return
    SparseAddress
    (
        std::move(rowPtr_),
        std::move(column_)
    );
}
}