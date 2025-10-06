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
    #include "createPhi.H"

    

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    #include "CourantNo.H"
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;
        dimensionedScalar dt= runTime.deltaT();                     // dt with unit	
        
        // phi is denfined on the faces of the mesh
	// fvc::interpolate(U) is Uf e.g linear
		    
	
	phi = fvc::interpolate(U) & mesh.Sf();
        fvScalarMatrix rhoEqn //Eqn 1 
	(
	   fvm::ddt(rho)
	+  fvc::div(phi,rho)
	);
	rhoEqn.solve(); 
        
        fvVectorMatrix rhoUEqn //Eqn 2 
	(
	   fvm::ddt(rhoU)
	+  fvc::div(phi,rhoU)
	==
	-  fvc::grad(p)
	);
	rhoUEqn.solve();
	U = rhoU / rho; //Eqn 3 
	
	
	fvScalarMatrix rhoEEqn //Eqn 4 
	(
	   fvm::ddt(rhoE)
	+  fvc::div(phi,rhoE)
	==
	-  fvc::div(phi,p)
	);
	rhoEEqn.solve();
	p = rho*R*T; //Eqn 5
	T = (rhoE - 0.5*rho*magSqr(U))/Cv/rho; //Eqn 6
	

        
        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
