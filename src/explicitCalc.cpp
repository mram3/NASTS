#include "Matrix/explicitCalc.h"

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

    const auto& faces = phi.getMesh().getFaces();
    const auto& cells = phi.getMesh().getCells();

    VectorField gradPhi("gradPhi", phi.getMesh());
    VectorField oldGradPhi("oldGradPhi", phi.getMesh());//to store previous corrected value

    //Correction Step is explicitly hardcoded to 2 following advice in the book
    //" Doing excessive number of iterations may cause oscillations and usually not more than 
    //  two iterations are performed."
    for(int correctionStep = 0; correctionStep < 3; ++correctionStep){

        for(const auto& cell: cells){
            gradPhi.vectorVariable[cell.cellID].x = 0.0;
            gradPhi.vectorVariable[cell.cellID].y = 0.0;
        }//initialized to zero so to dynamically add fluxes

        for(const auto& face : faces){
            int oID = face.ownerCellID;
            int nID = face.neighborCellID;

            //some Abbreviations: C     - primary cell's center
            //                    F     - neighbor cell's center
            //                    f     - face center
            //                    f'    - center of line joining two cell's center
            double xf = face.face_center_x;
            double yf = face.face_center_y;

            double Sx = face.face_normal_x*face.face_area;
            double Sy = face.face_normal_y*face.face_area;

            double xC = cells[oID].cell_center_x;
            double yC = cells[oID].cell_center_y;

            double xCf = xC - xf;// vector pointing from face to Cell center
            double yCf = yC - yf;

            if(face.physicalTag == -1){
                
                double xF = cells[nID].cell_center_x;
                double yF = cells[nID].cell_center_y;
                
                double xCF = xC - xF;// vector pointing from neighbor cell to primary cell
                double yCF = yC - yF;
                
                double q = - ((xCf*xCF) + (yCf*yCF)) / ((xCF*xCF) + (yCF*yCF));
                
                double xf_dash = xC + q * xCF;// locating f'
                double yf_dash = yC + q * yCF;

                double xFf_dash = xF - xf_dash;// vector from f' to neighbor cell
                double yFf_dash = yF - yf_dash;

                //geometric weighing factor for interpolation of value at f' from C and F
                double gC = (xFf_dash*xFf_dash + yFf_dash*yFf_dash) / (xCF*xCF + yCF*yCF);
                gC = std::sqrt(gC);

                double phiC = phi.scalarVariable[oID];
                double phiF = phi.scalarVariable[nID];
                double phif_dash = gC * phiC + (1 - gC) * phiF;// interpolating value

                //interpolating grad f' from grad at C and F for correcting gradient
                //by better estimating f from f'. 0 at first iteration, meaning f = f'.
                
                double xGradPhif_dash = gC * oldGradPhi.vectorVariable[oID].x 
                                        + (1 - gC) * oldGradPhi.vectorVariable[nID].x;
                double yGradPhif_dash = gC * oldGradPhi.vectorVariable[oID].y 
                                        + (1 - gC) * oldGradPhi.vectorVariable[nID].y;

                double xff_dash = xf - xf_dash;
                double yff_dash = yf - yf_dash;

                //Subsequent iteration we improve f = f' + grad f' * (rf - rf')
                double phif = phif_dash + xGradPhif_dash * xff_dash + yGradPhif_dash * yff_dash;

                gradPhi.vectorVariable[oID].x += phif*Sx;
                gradPhi.vectorVariable[oID].y += phif*Sy;
                gradPhi.vectorVariable[nID].x -= phif*Sx;
                gradPhi.vectorVariable[nID].y -= phif*Sy;
            }

            else{
                int tag = face.physicalTag;
                auto id = phi.boundaryConditionMap.find(tag);
                auto bcType = id->second.bc;
                double bcValue = id->second.value;

                if(bcType == BoundaryConditionType::Dirichlet){
                    double phif = bcValue;
                    //double phif = 2.0 * face.face_center_x + 3.0 * face.face_center_y; 
                    gradPhi.vectorVariable[oID].x += phif*Sx;
                    gradPhi.vectorVariable[oID].y += phif*Sy;
                }
                else if(bcType == BoundaryConditionType::Neumann){
                    double phiC = phi.scalarVariable[oID];

                    double normDist = std::abs(xCf*face.face_normal_x 
                                             + yCf*face.face_normal_y);

                    double phif = phiC + normDist*bcValue;

                    gradPhi.vectorVariable[oID].x += phif*Sx;
                    gradPhi.vectorVariable[oID].y += phif*Sy;
                }
            }
        }
        
        for(const auto& cell : cells){
            double volume = cell.cell_volume;

            gradPhi.vectorVariable[cell.cellID].x /= volume;
            gradPhi.vectorVariable[cell.cellID].y /= volume;

            oldGradPhi.vectorVariable[cell.cellID] = gradPhi.vectorVariable[cell.cellID];
        }
    }

    return gradPhi;
}

}