#include "Matrix/FvMatrix.h"

#include <algorithm>

void FvMatrix::addDiagCoeffs
(
    int cellID,
    double value
)
{
    int dIdx = A_.rowPtr()[cellID];//diagonals are stored at rowPtr[i] itself
                                   // as outlined in CsrAddress

    A_.values()[dIdx] += value;
}

void FvMatrix::addOffDiagCoeffs
(
    int faceID,
    double oFlux,
    double nFlux
)
{
    int oPtr = sparsityPattern_.oValPtr()[faceID];//the value pointer arrays stored in csrAddress
    int nPtr = sparsityPattern_.nValPtr()[faceID];

    A_.values()[oPtr] += oFlux;
    A_.values()[nPtr] += nFlux;
}

void FvMatrix::addSource
(
    int row,
    double value
)
{
    source_[row] += value;
}