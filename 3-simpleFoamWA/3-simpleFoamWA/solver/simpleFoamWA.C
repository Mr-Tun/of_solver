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
	#include "fvcGrad.H"          // 显式梯度计算     
	#include "bound.H" 

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
		// phi is denfined on the faces of the mesh
		// fvc::interpolate(U) is Uf e.g linear
                phi = fvc::interpolate(U) & mesh.Sf();                       //defined the convective flux		
                dimensionedScalar dt= runTime.deltaT();                     // dt with unit	
		volScalarField G = nut * dev(twoSymm(fvc::grad(U))) && T(fvc::grad(U));
		
		fvScalarMatrix epsilonEqn //Eqn 5
	        (
                        fvm::div(phi,epsilon)
                     -  fvm::laplacian(nut/sigmaEpsilon+nu, epsilon)
		     ==
		        c1 * G * epsilon / k
		     -  fvm::Sp(c2 * epsilon / k, epsilon)   
	        );  	
		epsilonEqn.relax();
		epsilonEqn.solve();
		bound(epsilon,dimensionedScalar(epsilon.dimensions(),0));    
		     
	        fvScalarMatrix kEqn //Eqn 6
	        (
	            
		        fvm::div(phi,k)
                     -  fvm::laplacian(nut/sigmaK+nu, k)
		     ==
		        G
		     -  fvm::Sp(epsilon/k,k)   
	        );  	
		    kEqn.relax();
		    kEqn.solve();
		    
		    nut = cu * magSqr(k) / epsilon;  //Eqn 7
		
		
		
	            
		    
		    fvVectorMatrix UEqn //Eqn 1 
		    (
			    fvm::div(phi,U)
			 -  fvm::laplacian(nu+nut, U)
			 ==
			 -  fvc::grad(p)
		    );
		    UEqn.relax();
		    UEqn.solve();
		    
		    volVectorField HbyA(constrainHbyA(UEqn.H()/UEqn.A(),U,p));  //Eqn 2
		    //volVectorField HbyA = UEqn.H()/UEqn.A();
		    
		    fvScalarMatrix pEqn //Eqn 3 
		    (
			 
			    fvm::laplacian(1/UEqn.A(),p)
			 ==
			    fvc::div(HbyA)
		    );
		    
		    pEqn.solve();
		    p.relax();
		    
		    U = U - 1/UEqn.A()*fvc::grad(p);  //Eqn 4
		    
		    
		    
		    
		}

		runTime.write();
	    }

	    Info<< "End\n" << endl;

	    return 0;
	}


	// ************************************************************************* //
