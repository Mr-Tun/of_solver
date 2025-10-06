volScalarField rhoOld = rho.oldTime();
volVectorField rhoUOld = rhoU.oldTime();
volScalarField rhoEOld = rhoE.oldTime();

rho = 0.75*rhoOld 
  + 0.25*(rho + runTime.deltaT()*(-fvc::div(phi, rho)));

rhoU = 0.75*rhoUOld 
  + 0.25*(rhoU + runTime.deltaT()*(-fvc::div(phi, rhoU) - fvc::grad(p)));

rhoE = 0.75*rhoEOld 
  + 0.25
   *(rhoE + runTime.deltaT()*(-fvc::div(phi, rhoE) - fvc::div(phi, p)));

T = (rhoE - 0.5*rho*magSqr(rhoU/rho))/Cv/rho;
p = rho*R*T;
U = rhoU/rho;
phi = fvc::interpolate(U) & mesh.Sf();

rho = 1.0/3.0*rhoOld 
  + 2.0/3.0*(rho + runTime.deltaT()*(-fvc::div(phi, rho)));

rhoU = 1.0/3.0*rhoUOld 
  + 2.0/3.0*(rhoU + runTime.deltaT()*(-fvc::div(phi, rhoU) - fvc::grad(p)));

rhoE = 1.0/3.0*rhoEOld 
  + 2.0/3.0
   *(rhoE + runTime.deltaT()*(-fvc::div(phi, rhoE) - fvc::div(phi, p)));

T = (rhoE - 0.5*rho*magSqr(rhoU/rho))/Cv/rho;
p = rho*R*T;
U = rhoU/rho;
phi = fvc::interpolate(U) & mesh.Sf();


