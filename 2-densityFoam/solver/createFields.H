IOdictionary physicalProperties
(
	IOobject
	(
	"physicalProperties",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
        IOobject::AUTO_WRITE
	)
);

dimensionedScalar R(physicalProperties.lookup("R"));
dimensionedScalar Cv = 5.0/2.0*R;



Info<< "Reading field of temperature\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field of pressure\n" << endl;
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

Info<< "Reading field of density\n" << endl;
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
    //mesh,
    //dimensionedScalar(dimless, 0.0)
    p/R/T,
    p.boundaryField().types()
);



Info<< "Reading field of velocity\n" << endl;
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

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    //mesh,
    //dimensionedScalar(dimless, 0.0)
    rho*U
);

Info<< "Reading field of density\n" << endl;
volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    //mesh,
    //dimensionedScalar(dimless, 0.0)
    rho*(Cv*T+0.5*magSqr(U))
);
