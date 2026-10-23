#include "Matrix/FvMatrix.h"
#include "Config.h"
#include "Algebra/Solvers/LinearSolvers.h"

#include <string>
#include <algorithm>
#include <memory>

//Constructors
FvMatrix::FvMatrix
(
    ScalarField& phi
)
:   CsrMatrix(phi.fvmesh().fvAddress()),
    phi_(phi),
    source_(phi.scalarVariable.size()),
    oValPtr_(phi.fvmesh().numInternalfaces(), -1),
    nValPtr_(phi.fvmesh().numInternalfaces(), -1)
{
    fillValPtrs();
}

FvMatrix::FvMatrix
(
    ScalarField& phi,
    std::shared_ptr<SparseAddress> address
)
:   CsrMatrix(address),
    phi_(phi),
    source_(phi.scalarVariable.size())
{

}

void FvMatrix::fillValPtrs()
{
    const auto& faces = phi_.fvmesh().faces();
    size_t numInFace = phi_.fvmesh().numInternalfaces();

    for(size_t i = 0; i < numInFace; ++i){
        int nID = faces[i].neighborCellID;
        int oID = faces[i].ownerCellID;

        const auto& col = sparsity_->column();
        const auto& rPt = sparsity_->rowPtr();

        //for each face we to owner cell row search where the neighbor cell's id is
        auto start = col.begin() + rPt[oID] + 1;
        auto end = col.begin() + rPt[oID + 1];

        auto valPtr = std::lower_bound(start, end, nID);

        oValPtr_[i] = std::distance(col.begin(), valPtr);

        //do the same for owner cell ID
        start = col.begin() + rPt[nID] + 1;
        end = col.begin() + rPt[nID + 1];

        valPtr = std::lower_bound(start, end, oID);

        nValPtr_[i] = std::distance(col.begin(), valPtr);
    }
}

void FvMatrix::addDiagCoeffs
(
    int cellID,
    double value
)
{
    int diagIdx = sparsity_->rowPtr()[cellID];//diagonals are stores at rowPtr[i] itself
    values_[diagIdx] += value;                 //as outlined in fvAddressing
}

void FvMatrix::addOffDiagCoeffs
(
    int faceID,
    double oVal,
    double nVal
)
{
    int oPtr = oValPtr_[faceID];
    int nPtr = nValPtr_[faceID];

    values_[oPtr] += oVal;
    values_[nPtr] += nVal;
}

void FvMatrix::addSource
(
    int cellID,
    double val
)
{
    source_[cellID] += val;
}

//The prime solve method of fv matrix
void FvMatrix::solve
(

)
{
    //name of the field eqn is solved for
    std::string fieldName = phi_.fieldName();

    //Choosing the solver for the specific field
    auto settings = Config::access().data()["LinearSolvers"][fieldName];

    //calling the solver
    auto solver = LinearSolvers::make(settings);

    solver->solve(*this, source_, phi_.scalarVariable);

}

//Global Solve Methods

void solve
(
    FvMatrix&& eqn
)
{
    eqn.solve();
}

//Global operators

FvMatrix& operator==
(
    FvMatrix& eqn,
    double rhs
)
{
    if(rhs == 0){
        return eqn;
    }
    
    for(size_t i = 0; i < eqn.source().size(); ++i){
        eqn.source()[i] += rhs;
    }
    
    return eqn;
}

FvMatrix&& operator==
(
    FvMatrix&& eqn,
    double rhs
)
{
    if(rhs == 0){
        return std::move(eqn);
    }
    
    for(size_t i = 0; i < eqn.source().size(); ++i){
        eqn.source()[i] += rhs;
    }
    
    return std::move(eqn);
}

FvMatrix operator-
(
    FvMatrix&& matA,
    FvMatrix&& matB
)
{
    if(matA.isDiag){
        const auto& rowPtrB = matB.pattern()->rowPtr();
        const auto& rowPtrA = matA.pattern()->rowPtr();
        const auto& valuesA = matA.values();
        auto& valuesB = matB.values();
        const auto& sourceA = matA.source();
        auto& sourceB = matB.source();

        for(size_t i = 0; i < matA.source().size(); ++i){
            sourceB[i] += sourceA[i];

            valuesB[rowPtrB[i]] += valuesA[rowPtrA[i]];
        }

        return std::move(matB);
    }

    else if(matB.isDiag){
        const auto& rowPtrA = matA.pattern()->rowPtr();
        const auto& rowPtrB = matB.pattern()->rowPtr();
        const auto& valuesB = matB.values();
        auto& valuesA = matA.values();
        const auto& sourceB = matB.source();
        auto& sourceA = matA.source();

        for(size_t i = 0; i < matB.source().size(); ++i){
            sourceA[i] += sourceB[i];

            valuesA[rowPtrB[i]] += valuesB[rowPtrA[i]];
        }

        return std::move(matA);
    }
}