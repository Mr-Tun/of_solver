IOdictionary physicalProperties
(
    IOobject 
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
Info<< "Reading constants \n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("d",dimDensity,physicalProperties.lookup("rho"))
);
volScalarField rhop
(
    IOobject
    (
        "rhop",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("d",dimDensity,physicalProperties.lookup("rhop"))
);
dimensionedScalar nu("d",dimensionSet(0, 2, -1, 0, 0, 0, 0),physicalProperties.lookup("nu"));
dimensionedScalar dp("d",dimLength,physicalProperties.lookup("dp"));
dimensionedScalar g("d",dimAcceleration,physicalProperties.lookup("g"));
dimensionedVector e("d",dimless,physicalProperties.lookup("e"));
volScalarField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    dp*dp*(rhop-rho)*g/18/nu/rho
);






Info<< "Reading field p\n" << endl;
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

Info<< "Reading field c\n" << endl;
volScalarField c
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
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

Info<< "Reading field Up\n" << endl;
volVectorField Up
(
    IOobject
    (
        "Up",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U + Us * e
);
Info<< "Reading field F\n" << endl;
volVectorField F
(
    IOobject
    (
        "F",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -3*3.1415926*nu*rho*dp*(U-Up)
);

#include "createPhi.H"


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
