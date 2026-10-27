#include "Applications/UnsteadyLaplacian.h"
#include "Matrix/FvMatrix.h"
#include "FiniteVolumeMethods/implicitAssembly/implicitAssembly.h"
#include "Controls/Time.h"
#include "Controls/SteadyControl.h"
#include "Config.h"
#include "Type.h"

#include <iostream>
#include <stdexcept>

UnsteadyLaplacian::UnsteadyLaplacian
(

)
: fvmesh_(),
  phi_(Config::access().data()["Fields"].begin().key(), fvmesh_)
{
    check();
    gamma_ = Config::access().data()["Gamma"];
    rho_   = Config::access().data()["Rho"];
}

void UnsteadyLaplacian::run
(

)
{
    createField();
    Time time;
    SteadyControl steadystate;

    double source = Config::access().data().value("Source", 0.0);

    int nNOC = Config::access().data().value("nNonOrthogonalCorrections", 0);

    while(time.isRunning()){

        for(int i = 0; i <= nNOC; ++i){

            solve(
                implicitAssembly::Ddt(rho_, phi_) - 
                implicitAssembly::laplacian(gamma_, phi_) == source
            );
        }

        std::cout << "Time : " << time.currTime() <<" | Temperature: " 
                  << phi_.scalarVariable[12] << "\n\n";

        if(steadystate.isReached(phi_))
        {
            std::cout << "Steady State is reached at time: " << time.currTime() << "\n";
            break;
        } 

        phi_.oldScalarVariable = phi_.scalarVariable;

    }
}

void UnsteadyLaplacian::createField
(

)
{

    auto fieldName = phi_.fieldName();
    auto fieldConfig = Config::access().data()["Fields"][fieldName];

    phi_.create(fieldConfig);
}