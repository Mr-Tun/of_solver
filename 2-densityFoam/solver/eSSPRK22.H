volScalarField rhoOld = rho.oldTime();
volVectorField rhoUOld = rhoU.oldTime();
volScalarField rhoEOld = rhoE.oldTime();

rho = 0.5*rhoOld 
  + 0.5*(rho + runTime.deltaT()*(-fvc::div(phi, rho)));

rhoU = 0.5*rhoUOld 
  + 0.5*(rhoU + runTime.deltaT()*(-fvc::div(phi, rhoU) - fvc::grad(p)));

rhoE = 0.5*rhoEOld 
  + 0.5
   *(rhoE + runTime.deltaT()*(-fvc::div(phi, rhoE) - fvc::div(phi, p)));

T = (rhoE - 0.5*rho*magSqr(rhoU/rho))/Cv/rho;
p = rho*R*T;
U = rhoU/rho;
phi = fvc::interpolate(U) & mesh.Sf();
        

