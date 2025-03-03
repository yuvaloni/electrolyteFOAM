Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
    IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    dimensionedScalar I
    (
        "I",
        dimensionSet(0, -2, 0, 0, 0, 0, 0),
        transportProperties
    );
    dimensionedScalar z1
    (
        "z1",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        transportProperties
    );

    dimensionedScalar z2
    (
        "z2",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        transportProperties
    );

    dimensionedScalar zeta1
    (
        "zeta1",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        transportProperties
    );

    dimensionedScalar zeta2
    (
        "zeta2",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        transportProperties
    );

    dimensionedScalar e
    (
        "e",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        transportProperties
    );

    dimensionedScalar kappa_d
    (
        "kappa_d",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        transportProperties
    );

    dimensionedScalar eta
    (
        "eta",
        dimensionSet(0, 2, 0, 0, 0, 0, 0),
        transportProperties
    );

    dimensionedScalar beta
    (
        "beta",
        dimensionSet(0, -2, 1, 0, 0, 0, 0),
        transportProperties
    );

    dimensionedScalar rho
    (
        "rho",
        dimensionSet(0, 2, 0, 0, 0, 0, 0),
        transportProperties
    );

    dimensionedScalar D1
    (
        "D1",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        transportProperties
    );

    dimensionedScalar D2
    (
        "D2",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        transportProperties
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





Info<< "Reading field psi\n" << endl;
volScalarField psi
(
    IOobject
    (
        "Psi",
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

Info<< "Reading field n1\n" << endl;
volScalarField n1
(
    IOobject
    (
        "n1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field n2\n" << endl;
volScalarField n2
(
    IOobject
    (
        "n2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field v1\n" << endl;
volVectorField v1
(
    IOobject
    (
        "v1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field v2\n" << endl;
volVectorField v2
(
    IOobject
    (
        "v2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

#include "createFvOptions.H"
