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




Info<< "Reading field k\n" << endl;
volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field epsilon\n" << endl;
volScalarField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading physicalProperties\n" << endl;

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
dimensionedScalar c1 = 1.44;
dimensionedScalar c2 = 1.92;
dimensionedScalar cu = 0.09;
dimensionedScalar sigmaK = 1.00;
dimensionedScalar sigmaEpsilon = 1.30;

Info<< "Reading field nut\n" << endl;
volScalarField nut
(
    IOobject
    (
        "nut",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    //mesh
     cu * k * k / epsilon
);

Info<< "Reading diffusivity nu\n" << endl;

dimensionedScalar nu // viscousity
(
    physicalProperties.lookup("nu")
);

#include "createPhi.H"

#include "createFvModels.H"
#include "createFvConstraints.H"
