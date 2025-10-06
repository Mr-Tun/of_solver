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
    pisoFoam

Description
    Transient solver for incompressible, turbulent flow, using the PISO
    algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "viscosityModel.H"
#include "incompressibleMomentumTransportModels.H"
#include "pisoControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        #include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            fvModels.correct();
            fvVectorMatrix UEqn
            (
               fvm::ddt(rho,U)
             + fvm::div(phi,U) //U_f & mesh.Sf() 不可压缩的phi，现在算可压缩要加上密度
             - fvm::laplacian(mu,U)   
            ==
             -fvc::grad(p)             
            );
            UEqn.relax();
            UEqn.solve();
            
            
            //#include "UEqn.H"//U*=U ^t
             
            // --- PISO loop
            while (piso.correct())
            {
                  
                //#include "pEqn.H"
                volScalarField A (UEqn.A());
                volVectorField HbyA(constrainHbyA(UEqn.H()/A, U, p)); 
                surfaceScalarField phip = fvc::interpolate(psi*HbyA) & mesh.Sf();
                fvScalarMatrix pEqn
                (
                   fvm::ddt(psi,p)
                +  fvc::div(HbyA*(rho0-psi*p0))
                +  fvc::div(phip,p)
                -  fvm::laplacian(rho/A,p)
                   
                );
                pEqn.relax();
                pEqn.solve();
                rho = rho0 + psi * (p+p0);
            }
        }

        viscosity->correct();
        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
