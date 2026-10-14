#include "Matrix/CsrAddress.h"

#include <algorithm>

CsrAddress::CsrAddress
(
    const Grid& mesh
)
{
    auto const& cells = mesh.getCells();
    auto const& faces = mesh.getFaces();

    int nCells = cells.size();
    int numInFace = mesh.internalFaces();

    //rowPtr is nothing but a count of non-zero elements before that row
    rowPtr_.assign(nCells + 1, 0);

    //so we lookup for number of non-zero contribution by each cell we have, by looking at how
    //many internal faces are shared by each cell and ofcourse a cell will have it's own
    //non-zero coefficient aP...

    std::vector<int> nnz(nCells, 1); //...so we start with 1

    for(int i = 0; i < numInFace; ++i){//boundary faces do not contribute to matrix A
        auto& face = faces[i];

        nnz[face.ownerCellID]++; 
        nnz[face.neighborCellID]++; 
    }

    //now that we have for each cell how many nnz it has, now we just progressively add it
    //to get nnz before each row

    rowPtr_[0] = 0; //has to be 
    for(int i = 0; i < nCells; ++i){
        rowPtr_[i+1] = rowPtr_[i] + nnz[i];
    }

    column_.assign(rowPtr_[nCells], 0);

    std::vector<int> tracker = rowPtr_; //to track within a row
    for(int i = 0; i < nCells; ++i){//just the cell's own column
        column_[tracker[i]] = cells[i].cellID;
        tracker[i]++;
    }

    for(int i = 0; i < numInFace; ++i){
        int nID = faces[i].neighborCellID;
        int oID = faces[i].ownerCellID;

        column_[tracker[oID]] = nID;//each face's owner cell will have neighbor ID column
        column_[tracker[nID]] = oID;//and vice versa

        tracker[oID]++;
        tracker[nID]++;
    }

    for(int i = 0; i < nCells; ++i){//keeping the diagonal at first position and sorting others
        std::sort
        (
            column_.begin() + rowPtr_[i]+1, column_.begin() + rowPtr_[i+1]
        );
    }

    //pointer arrays to values array of CSR Matrix
    oValPtr_.assign(numInFace, 0);
    nValPtr_.assign(numInFace, 0);

    for(int i = 0; i < numInFace; ++i){
        int nID = faces[i].neighborCellID;
        int oID = faces[i].ownerCellID;

        //for each face we go to owner row, see where neighbor column is, store it
        auto start = column_.begin() + rowPtr_[oID] + 1;
        auto end   = column_.begin() + rowPtr_[oID + 1];

        auto valPtr = std::lower_bound(start, end, nID);

        oValPtr_[i] = std::distance(column_.begin(), valPtr);

        //same for neighbor row of the face
        start = column_.begin() + rowPtr_[nID] + 1;
        end = column_.begin() + rowPtr_[nID + 1];

        valPtr = std::lower_bound(start, end, oID);

        nValPtr_[i] = std::distance(column_.begin(), valPtr);
    }
}