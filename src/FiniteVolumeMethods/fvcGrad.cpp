#include "FiniteVolumeMethods/explicitCalc.h"

#include <cmath>

namespace explicitCalc
{

VectorField grad
(
    const ScalarField& phi
)
{
    //brief Explicit Calculation of Gradient at Cell centers based on previous iteration values
    //      Ref : The FVM in CFD by F. Moukalled, L. Mangani, M. Darwish
    //      Gradient calculation using Green-Gauss method with skewness correction
    //      For skewness correction, method 3 in the book is adopted

    const auto& faces = phi.fvmesh().faces();
    const auto& cells = phi.fvmesh().cells();
    int numInFace = phi.fvmesh().numInternalfaces();
    int numFaces  = faces.size();

    VectorField gradPhi("gradPhi", phi.fvmesh());
    //VectorField oldGradPhi("oldGradPhi", phi.fvmesh());//to store previous corrected value

    //Correction Step is explicitly hardcoded to 2 following advice in the book
    //" Doing excessive number of iterations may cause oscillations and usually not more than 
    //  two iterations are performed."
    for(int correctionStep = 0; correctionStep < 3; ++correctionStep){

        for(const auto& cell: cells){
            gradPhi.vectorVariable[cell.cellID].x = 0.0;
            gradPhi.vectorVariable[cell.cellID].y = 0.0;
        }//initialized to zero so to dynamically add fluxes

        for(int i = 0; i < numInFace; ++i){
            const auto& face = faces[i];
            int oID = face.ownerCellID;
            int nID = face.neighborCellID;

            //some Abbreviations: C     - primary cell's center
            //                    F     - neighbor cell's center
            //                    f     - face center
            //                    f'    - center of line joining two cell's center
            //To match with varaible names used in the book

            double xCf = face.xOf;// vector pointing from face f to Cell center C
            double yCf = face.yOf;

            double xCF = face.xON;// vector pointing from neighbor cell to primary cell C
            double yCF = face.yON;

            double q = - ((xCf*xCF) + (yCf*yCF)) / ((xCF*xCF) + (yCF*yCF));

            //geometric weighing factor for interpolation of value at f' from C and F
            double gC = 1 + q;

            double phiC = phi.scalarVariable[oID];
            double phiF = phi.scalarVariable[nID];
            double phif_dash = gC * phiC + (1 - gC) * phiF;// interpolating value

            //interpolating grad f' from grad at C and F for correcting gradient
            //by better estimating f from f'. 0 at first iteration, meaning f = f'.
                
            double xGradPhif_dash = gC * gradPhi.oldVectorVariable[oID].x 
                                    + (1 - gC) * gradPhi.oldVectorVariable[nID].x;
            double yGradPhif_dash = gC * gradPhi.oldVectorVariable[oID].y 
                                    + (1 - gC) * gradPhi.oldVectorVariable[nID].y;

            double xff_dash = xCf + q * xCF;
            double yff_dash = yCf + q * yCF;

            //Subsequent iteration we improve f = f' + grad f' * (rf - rf')
            double phif = phif_dash + xGradPhif_dash * xff_dash + yGradPhif_dash * yff_dash;

            gradPhi.vectorVariable[oID].x += phif*face.face_area_x;
            gradPhi.vectorVariable[oID].y += phif*face.face_area_y;
            gradPhi.vectorVariable[nID].x -= phif*face.face_area_x;
            gradPhi.vectorVariable[nID].y -= phif*face.face_area_y;
        }

        for(int b = numInFace; b < numFaces; ++b){
            const auto& face = faces[b];
            int tag = face.physicalTag;
            auto id = phi.boundaryConditionMap.find(tag);
            auto bcType = id->second.bc;
            double bcValue = id->second.value;

            int oID = face.ownerCellID;

            if(bcType == BoundaryConditionType::Dirichlet){
                double phif = bcValue;
                //double phif = 2.0 * face.face_center_x + 3.0 * face.face_center_y; 
                gradPhi.vectorVariable[oID].x += phif*face.face_area_x;
                gradPhi.vectorVariable[oID].y += phif*face.face_area_y;
            }
            else if(bcType == BoundaryConditionType::Neumann){
                double phiC = phi.scalarVariable[oID];

                double xCf = face.xOf;
                double yCf = face.yOf;
                double normDist = std::abs(xCf*face.face_normal_x + yCf*face.face_normal_y);

                double phif = phiC + normDist*bcValue;

                gradPhi.vectorVariable[oID].x += phif*face.face_area_x;
                gradPhi.vectorVariable[oID].y += phif*face.face_area_y;
            }
        }

        for(const auto& cell : cells){
            double volume = cell.cell_volume;

            gradPhi.vectorVariable[cell.cellID].x /= volume;
            gradPhi.vectorVariable[cell.cellID].y /= volume;

            gradPhi.oldVectorVariable[cell.cellID] = gradPhi.vectorVariable[cell.cellID];
        }
    }
        
    return gradPhi;
}
  
}