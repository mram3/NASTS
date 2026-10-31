#include "Applications/ScalarTransport.h"
#include "Matrix/FvMatrix.h"
#include "FiniteVolumeMethods/implicitAssembly/implicitAssembly.h"
#include "Controls/Time.h"
#include "Controls/SteadyControl.h"
#include "Controls/Output.h"
#include "Config.h"
#include "Type.h"

#include <iostream>
#include <stdexcept>

ScalarTransport::ScalarTransport
(

)
: fvmesh_(),
  time_(),
  output_(fvmesh_, time_)
{
    check();
    gamma_ = Config::access().data()["Gamma"];
    rho_   = Config::access().data()["Rho"];

    std::string fieldName = "T";//Config::access().data()["Fields"]["T"];
    phi_ = std::make_unique<ScalarField>(fieldName, fvmesh_);

    std::string velName = "U";
    U_ = std::make_unique<VectorField>(velName, fvmesh_);

    std::string fluxName = "Flux";
    flux_ = std::make_unique<SurfaceScalarField>(fluxName, fvmesh_);
}

void ScalarTransport::run
(

)
{
    createField();
    SteadyControl steadystate;

    double source = Config::access().data().value("Source", 0.0);

    int nNOC = Config::access().data().value("nNonOrthogonalCorrections", 0);

    //storing the initial condition
    output_.writeVtk(true);

    while (time_.isRunning())
    {

        phi_->storeField();

        for(int i = 0; i <= nNOC; ++i){

            solve(
                implicitAssembly::Ddt(rho_, *phi_)
              + implicitAssembly::div(*flux_, *phi_)
              - implicitAssembly::laplacian(gamma_, *phi_) == source
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

void ScalarTransport::createField
(

)
{
    auto fieldName = phi_->fieldName();
    auto fieldConfig = Config::access().data()["Fields"][fieldName];

    phi_->create(fieldConfig);

    auto velConfig = Config::access().data()["Fields"]["U"];

    U_->createVelocityField(velConfig);

    flux_->calcFlux(rho_, *U_);

    output_.registerScalar(*phi_);
    output_.registerVector(*U_);
}