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
	#include "fvModels.H"
	#include "fvConstraints.H"
	#include "simpleControl.H"

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

	    #include "CourantNo.H"

	    while (simple.loop(runTime))
	    {
		Info<< "Time = " << runTime.userTimeName() << nl << endl;

		fvModels.correct();

		while (simple.correctNonOrthogonal())
		{
                    dimensionedScalar dt= runTime.deltaT();                     // dt with unit	
		    fvVectorMatrix UEqn //Eqn 1 
		    (
			    fvm::ddt(U)
			 +  fvc::div(phi,U)
			 -  fvc::laplacian(DT, U)
		    );
		    UEqn.solve();
		    
		    // phi is denfined on the faces of the mesh
		    // fvc::interpolate(U) is Uf e.g linear
		    
		     phi = fvc::interpolate(U) & mesh.Sf();
		    
		    fvScalarMatrix pEqn //Eqn 2 
		    (
		       
		        fvm::laplacian(dt,p) == fvc::div(phi)  
		    );
		    
		    pEqn.solve();
		    
		     //Eqn 3
		   
		     U = -dt * fvc::grad(p) + U; 
		    
		}

		runTime.write();
	    }

	    Info<< "End\n" << endl;

	    return 0;
	}


	// ************************************************************************* //
