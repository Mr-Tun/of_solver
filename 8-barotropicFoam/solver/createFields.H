//参考物理量
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
   mesh
);
dimensionedScalar  rho0(dimDensity,1);
//dimensionedScalar  p0(dimPressure,100000);
dimensionedScalar  p0(dimensionSet(1,-1,-2,0,0),100000);
dimensionedScalar  mu(rho0.dimensions()*dimViscosity,1e-5);
Info<< "Reading field p\n" << endl;



volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
   mesh, 
   dimensionedScalar(dimensionSet(0,-2,2,0,0),sqr(1.0/340.0))
);



volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho0 + psi*(p - p0) //psi 可压缩性
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//#include "createPhi.H"
//surfaceScalarField phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf() );
surfaceScalarField phi ("phi",fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf() ));
pressureReference pressureReference(p, piso.dict());

mesh.schemes().setFluxRequired(p.name());


autoPtr<viscosityModel> viscosity(viscosityModel::New(mesh));

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, viscosity)
);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
