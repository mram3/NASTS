#include "Applications/Laplacian.h"
#include "Matrix/FvMatrix.h"
#include "FiniteVolumeMethods/implicitAssembly.h"
#include "Config.h"
#include "Type.h"

#include <iostream>
#include <stdexcept>

Laplacian::Laplacian
(

) 
:   mesh_(),
    phi_(Config::access().data()["Fields"].begin().key(), mesh_)
{
    gamma_ = Config::access().data()["Gamma"];
}

void Laplacian::run
(

)
{
    createField();
    
    double source = Config::access().data().value("Source", 0.0);

    int nNOC = Config::access().data().value("nNonOrthogonalCorrections", 0);

    for(int i = 0; i <= nNOC; ++i){

        FvMatrix TEqn = 
        (
            implicitAssembly::laplacian(gamma_, phi_)
        );

        TEqn == source;

        TEqn.solve();
    }

    std::cout << "Solution successfully converged\n";
    std::cout << "Temperature at center : " << phi_.scalarVariable[12] << "\n";
}

void Laplacian::createField
(

)
{

    auto fieldName = phi_.fieldName();
    auto fieldConfig = Config::access().data()["Fields"][fieldName];

    phi_.create(fieldConfig);
}