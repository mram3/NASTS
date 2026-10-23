#pragma once

#include "Matrix/SparseAddress.h"
#include "Grid/FvMesh.h"

namespace FvAddressing
{
    //From FvMesh, building the sparsity pattern
    SparseAddress buildFvAddress(const FvMesh& mesh);

    SparseAddress buildFvDiagAddress(const FvMesh& mesh);
}