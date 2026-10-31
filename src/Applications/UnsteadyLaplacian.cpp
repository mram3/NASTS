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
  time_(),
  output_(fvmesh_, time_)
{
    check();
    gamma_ = Config::access().data()["Gamma"];
    rho_   = Config::access().data()["Rho"];

    std::string fieldName = Config::access().data()["Fields"].begin().key();
    phi_ = std::make_unique<ScalarField>(fieldName, fvmesh_);
}

void UnsteadyLaplacian::run
(

)
{
    createField();
    SteadyControl steadystate;

    double source = Config::access().data().value("Source", 0.0);

    int nNOC = Config::access().data().value("nNonOrthogonalCorrections", 0);

    //storing the initial condition
    output_.writeVtk(true);

    while(time_.isRunning()){

        phi_->storeField();

        for(int i = 0; i <= nNOC; ++i){

            solve(
                implicitAssembly::Ddt(rho_, *phi_) - 
                implicitAssembly::laplacian(gamma_, *phi_) == source
            );
        }

        if(steadystate.isReached(*phi_))
        {
            if(time_.isSteadySolve()){
                std::cout << "Steady State is reached at iteration: " << time_.currIter() << "\n";
            }
            else{
                std::cout << "Steady State is reached at time: " << time_.currTime() << "\n";
            }
            //storing the steady state convergence
            output_.writeVtk(true);
            break;
        }

        output_.writeVtk();
    }
}

void UnsteadyLaplacian::createField
(

)
{

    auto fieldName = phi_->fieldName();
    auto fieldConfig = Config::access().data()["Fields"][fieldName];

    phi_->create(fieldConfig);

    output_.registerScalar(*phi_);
}