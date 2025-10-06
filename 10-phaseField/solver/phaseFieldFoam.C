/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    scalarTransportFoam

Description
    Solves the steady or transient transport equation for a passive scalar.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "fvcGrad.H"          // 显式梯度计算   
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;
    surfaceScalarField phi = fvc::interpolate(U) & mesh.Sf();                       //defined the convective flux	
    #include "CourantNo.H"
   
    Info<< "continuity error: " << max(fvc::div(phi)) << nl;
    //Sf magSf N
    const surfaceVectorField Sf = mesh.Sf();
    const surfaceScalarField magSf = mesh.magSf();             
    const surfaceVectorField N = Sf/magSf;
    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;
        /*fvScalarMatrix TEqn
        (
        fvm::ddt(T)
        +fvm::div(phi,T)
        );
        TEqn.relax();
        TEqn.solve();
        */
        volVectorField gradT = fvc::grad(T);
        volVectorField n = gradT/ (mag(gradT) + deltaN);
        surfaceVectorField nf = fvc::interpolate(n);
        surfaceScalarField negTf = fvc::interpolate(1-T);
        surfaceScalarField Gammaf = fvc::interpolate(mag(U));
        surfaceScalarField phiT = negTf*Gammaf*(nf&mesh.Sf());
        fvScalarMatrix TEqn
        (
          fvm::ddt(T)
        //+ U&fvm::grad(T)
        + fvm::div(phi,T)            
        - fvc::div(phi)*T 
        + fvm::div(phiT,T)  
        
        );
        TEqn.relax();
        TEqn.solve();
        
        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
