#include "Matrix/implicitAssembly.h"
#include "Matrix/explicitCalc.h"
#include "Fields/ScalarField.h"
#include "Fields/VectorField.h"

#include <stdexcept>
#include <string>
#include <cmath>

namespace implicitAssembly
{

FvMatrix laplacian
(
    double gamma,
    ScalarField& phi
)
{
    //initializing the matrix arrays
    FvMatrix Mat(phi);

    auto const& faces = phi.getMesh().getFaces();
    auto const& cells = phi.getMesh().getCells();

    //brief Implementation of diffusion on Unstructured Meshes
    //      Ref :- NPTEL: CFD using FVM by Kameswararao Anupindi
    //      Flux leaving cell 0 to cell 1 through face 1 (let's say) thus if we evaluate flux
    //      at faces it becomes coefficients for both cells. This suggests associating coeffs
    //      with cells rather than the faces

    VectorField gradPhi = explicitCalc::grad(phi);// for S.G contribution
    for(auto const& face : faces){
        if(face.physicalTag == -1){//looping through the internal faces

            int oID = face.ownerCellID;
            int nID = face.neighborCellID;

            double xC = cells[oID].cell_center_x;
            double yC = cells[oID].cell_center_y;
            double xF = cells[nID].cell_center_x;
            double yF = cells[nID].cell_center_y;

            double dx = xF - xC;
            double dy = yF - yC;

            double denominator = (face.face_normal_x*dx + face.face_normal_y*dy);

            //brief "Over-Relaxed Approach"
            //      a_nb = (gamma_f/d_zeta)*(Af * Af)/(Af * e_zeta)
            //      after simplification, a_nb = gamma_f * Af / nf * e_zeta
            double a_nb = gamma*face.face_area/(std::abs(denominator));

            Mat.addCoeffs(oID, oID, a_nb);
            Mat.addCoeffs(oID, nID, -a_nb);
            Mat.addCoeffs(nID, nID, a_nb);
            Mat.addCoeffs(nID, oID, -a_nb);

            //Secondary Gradient
            // T = Af - E (first evaluate T)

            double xA = face.face_area * face.face_normal_x;
            double yA = face.face_area * face.face_normal_y;
            double xE = a_nb * dx / gamma;//E is just the geometric part of a_nb
            double yE = a_nb * dy / gamma;//just projecting it and dividing by gamma 
            
            double xT = xA - xE;
            double yT = yA - yE;

            // (dphi/deta)f = T * grad(phi)f
            // But grad(phi) is only available at cell centers thus have to interpolated to face

            //brief Interpolating Gradients to Faces (section 9.4)
            //      Ref : The FVM in CFD by F. Moukalled, L. Mangani, M. Darwish
            //      "It is important for the stencil of the gradient at the face to be heavily
            //       based on the nodes straddling the face, which is not guarenteed by simple
            //       averaging practice."

            //first calculate geometric weight factor
            double xf = face.face_center_x;
            double yf = face.face_center_y;
            double xFf = xF - xf;
            double yFf = yF - yf;

            double gC = ((xFf * xFf) + (yFf * yFf)) / ((dx * dx) + (dy * dy));
            gC = std::sqrt(gC);

            //interpolating gradient values to face (crude simple average estimate)
            double xGradPhif = gC * gradPhi.vectorVariable[oID].x 
                             + (1 - gC) * gradPhi.vectorVariable[nID].x;
                             
            double yGradPhif = gC * gradPhi.vectorVariable[oID].y 
                             + (1 - gC) * gradPhi.vectorVariable[nID].y;

            //correction term : [phiF - phiC]/dCF - [(gradPhif)crude * eCF] eCF

            double phiC = phi.scalarVariable[oID];
            double phiF = phi.scalarVariable[nID];
            double dCF = std::sqrt(dx*dx + dy*dy);

            double nDiff = (phiF - phiC)/dCF; //[phiF - phiC]/dCF

            double projCrossDiff = (xGradPhif * dx + yGradPhif * dy)/dCF; //(gradPhif)crude*eCF

            double xCorr = (nDiff - projCrossDiff) * dx/dCF;
            double yCorr = (nDiff - projCrossDiff) * dy/dCF;

            //now update the correction
            xGradPhif += xCorr;
            yGradPhif += yCorr;

            //now computing S.G as gamma * (dphi/deta)f
            double sG = gamma * (xT * xGradPhif + yT * yGradPhif);

            //add it to the source term (rhs)
            Mat.addSource(oID, sG);
            Mat.addSource(nID, -sG);
        }

        else{
            int oID = face.ownerCellID;
            int tag = face.physicalTag;

            auto id = phi.boundaryConditionMap.find(tag);

            if(id == phi.boundaryConditionMap.end())
            {
                std::string errMsg = "Error: Missing Boundary Condition for physical tag "
                                      + std::to_string(tag);
                
                throw std::runtime_error(errMsg);
            }

            auto bcType = id->second.bc;
            auto bcValue = id->second.value;

            double dx = face.face_center_x - cells[oID].cell_center_x;
            double dy = face.face_center_y - cells[oID].cell_center_y;

            double denominator = (face.face_normal_x*dx + face.face_normal_y*dy);

            if(bcType == BoundaryConditionType::Dirichlet){
                double a_nb = gamma*face.face_area/(std::abs(denominator));
                //bcValue = 2 * face.face_center_x + 3 * face.face_center_y;

                //boundary cross diffusion (to be done if boundary cells are skewed)
                // A = E + T (T is )
                double xA = face.face_area * face.face_normal_x;
                double yA = face.face_area * face.face_normal_y;
                double xE = a_nb * dx / gamma;
                double yE = a_nb * dy / gamma; 
    
                double xT = xA - xE;
                double yT = yA - yE;

                //unlike interior faces, boundary faces do not have neighbor cells
                //thus cannot be interpolated, we take owner cell's gradient
                double xGradPhif = gradPhi.vectorVariable[oID].x;
                double yGradPhif = gradPhi.vectorVariable[oID].y;

                //now computing S.G as gamma * (dphi/deta)f
                double sG = gamma * (xT * xGradPhif + yT * yGradPhif);

                //b = add to the source term
                double b = (a_nb * bcValue) + sG;

                Mat.addCoeffs(oID, oID, a_nb);
                Mat.addSource(oID, b);
            }
            /*if(bcType == BoundaryConditionType::Dirichlet){

                double a_nb = gamma*face.face_area/(std::abs(denominator));
                bcValue = 2 * face.face_center_x + 3 * face.face_center_y;
                double b = a_nb*bcValue;

                Mat.addCoeffs(oID, oID, a_nb);
                Mat.addSource(oID, b);
            }*/
            else if(bcType == BoundaryConditionType::Neumann){
                double b = gamma*face.face_area*bcValue;
                Mat.addSource(oID, b);
            }
        }
    }

    return Mat;
}

}