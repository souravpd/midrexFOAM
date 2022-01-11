Info<< "Reading velocity field USolid\n" << endl;
volVectorField USolid
(
    IOobject
    (
        "USolid",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

USolid = dimensionedVector(USolid.dimensions(), Zero);

surfaceScalarField phiSolid
(
    IOobject
    (
        "phiSolid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(USolid)
);

if (args.optionFound("initialiseUBCs"))
{
    USolid.correctBoundaryConditions();
    phiSolid = fvc::flux(USolid);
}
// Default name for the pressure field
word pName("pSolid");

// Update name of the pressure field from the command-line option
args.optionReadIfPresent("pName", pName);

// Infer the pressure BCs from the velocity
wordList pBCTypes
(
    USolid.boundaryField().size(),
    fixedValueFvPatchScalarField::typeName
);

forAll(USolid.boundaryField(), patchi)
{
    if (USolid.boundaryField()[patchi].fixesValue())
    {
        pBCTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
    }
}

Info<< "Constructing Solid pressure field " << pName << nl << endl;
volScalarField pSolid
(
    IOobject
    (
        pName,
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(pName, sqr(dimVelocity), 0),
    pBCTypes
);

// Infer the velocity potential BCs from the pressure
wordList PhiBCTypes
(
    pSolid.boundaryField().size(),
    zeroGradientFvPatchScalarField::typeName
);

forAll(pSolid.boundaryField(), patchi)
{
    if (pSolid.boundaryField()[patchi].fixesValue())
    {
        PhiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
}

Info<< "Constructing velocity potential field Phi\n" << endl;
volScalarField PhiSolid
(
    IOobject
    (
        "PhiSolid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength*dimVelocity, 0),
    PhiBCTypes
);

Info<< "Reading solid field TSolid\n" << endl;
volScalarField TSolid
(
    IOobject
    (
        "TSolid",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading diffusivity DTSolid\n" << endl;

dimensionedScalar DTSolid
(
    transportProperties.lookup("DTSolid")
);

dimensionedScalar heffSolid
(
    transportProperties.lookup("heffSolid")
);