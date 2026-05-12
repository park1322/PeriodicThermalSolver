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

Info<< "Reading field T\n" << endl;
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

// volScalarField theta
// (
//     IOobject
//     (
//         "theta",
//         runTime.timeName(),
//         mesh
//     ),
//     T - Twall
// );

volScalarField lambda
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    // dimensionedScalar("lambda", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.0),
    // T.boundaryField().types()
);

// volScalarField lambda_x
// (
//     IOobject
//     (
//         "lambda_x",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     sqr(lambda)
// );

volScalarField lambda_x
(
    IOobject
    (
        "lambda_x",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless/sqr(dimLength), 0.0),  // 1/L^2
    lambda.boundaryField().types()                           // <- inherit types
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// START custom code 

// thermal diffusivity
Info<< "Reading value Dt\n" << endl;
dimensionedScalar Dt("Dt", laminarTransport);
Info<< "Dt is " << Foam::name(Dt.value()) << endl;

// END custom code

#include "createMRF.H"
#include "createFvOptions.H"
