#include "Tests.h"
#include "Grid/Grid.h"
#include "Grid/MeshReader.h"
#include "Grid/FvMesh.h"
#include "FiniteVolumeFields/ScalarField.h"
#include "Matrix/SparseAddress.h"
#include "Matrix/FvAddressing.h"
#include "Matrix/CsrMatrix.h"
#include "Matrix/FvMatrix.h"

#include <iostream>
#include <vector>
#include <algorithm>

bool isDiagFirst(const std::vector<int>& rowPtr, const std::vector<int>& col);

bool isColSorted(const std::vector<int>& rowPtr, const std::vector<int>& col);

bool isValPtrValid
(
    int n, const std::vector<Face>& faces,
    const std::vector<int>& oValPtr, 
    const std::vector<int>& nValPtr,
    const std::vector<int>& col
);

int main(){
    std::cout << "--FV Addressing Test--\n";

    try
    {
        //setting up the address
        Grid mesh = MeshReader::readGmshAscii22("tests/assets/structured.msh");

        FvMesh fvmesh(mesh);

        SparseAddress address = FvAddressing::buildFvAddress(fvmesh);

        const auto& rowPtr = address.rowPtr();
        const auto& col = address.column();

        //tests for address
        EXPECT_EQ
        (
            col.size(),
            (fvmesh.cells().size() + 2*fvmesh.numInternalfaces()),
            "Column Size check"
        );

        EXPECT_TRUE
        (
            isDiagFirst(rowPtr, col),
            "Diagonals are first in each row"
        );

        EXPECT_TRUE
        (
            isColSorted(rowPtr, col),
            "Off-Diagonal columns are sorted"
        );

        EXPECT_EQ
        (
            rowPtr[1] - rowPtr[0],
            3,
            "Boundary Address is fine"
        );

        EXPECT_EQ
        (
            rowPtr[12] - rowPtr[11],
            5,
            "Interiror Address is fine"
        );

        //test CSR Matrix' SpMV
        CsrMatrix csrMat(std::make_shared<SparseAddress>(address));
        std::vector<double> x(rowPtr.size() - 1, 1.0);
        std::vector<double> result(rowPtr.size() - 1, 0.0);

        csrMat.values().assign(col.size(), 1);

        csrMat.SpMV(x, result);

        EXPECT_TRUE(
            result[0] == 3 && 
            result[1] == 4 &&
            result[12] == 5,
            "Sparse Matrix-Vector is fine"
        );

        //test FV Matrix' Value Pointer Arrays
        ScalarField testPhi("test", fvmesh);
        
        FvMatrix fvMat(testPhi);

        EXPECT_TRUE
        (
            isValPtrValid
            (
                fvmesh.numInternalfaces(),
                fvmesh.faces(),
                fvMat.oValPtr(),
                fvMat.nValPtr(),
                col
            ),
            "Value Pointers are fine"
        );
    }
    catch(const std::exception& e)
    {
        std::cout << "\033[31m[ CRASH ]\033[0m " << e.what() << "\n";
        return 1;
    }

    return Tests::summary();
}

bool isDiagFirst
(
    const std::vector<int>& rowPtr, 
    const std::vector<int>& col
)
{
    for(int i = 0; i < rowPtr.size() - 1; ++i){
        if(col[rowPtr[i]] != i) return false;
    }

    return true;
}

bool isColSorted
(
    const std::vector<int>& rowPtr, 
    const std::vector<int>& col
)
{
    for(int i = 0; i < rowPtr.size() - 1; ++i){
        if(!(std::is_sorted
        (
            col.begin() + rowPtr[i] + 1,
            col.begin() + rowPtr[i+1]
        )))
        {
            return false;
        }
    }

    return true;
}

bool isValPtrValid
(
    int n, const std::vector<Face>& faces,
    const std::vector<int>& oValPtr, 
    const std::vector<int>& nValPtr,
    const std::vector<int>& col
)
{
    for(int i = 0; i < n; ++i){
        if(
            col[oValPtr[i]] != faces[i].neighborCellID ||
            col[nValPtr[i]] != faces[i].ownerCellID
        ){
            return false;
        }
    }

    return true;
}