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

label PhiRefCell = 0;
scalar PhiRefValue = 0;
setRefCell
(
    PhiSolid,
    simple.dict(),
    PhiRefCell,
    PhiRefValue
);
mesh.setFluxRequired(PhiSolid.name());
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








//reading fluid fields
Info<< "Reading field pFluid\n" << endl;
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

Info<< "Reading field UFluid\n" << endl;
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

Info<< "calculating face flux field phiFluid\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

Info<< "Reading fluid field T\n" << endl;
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


label pFluidRefCell = 0;
scalar pFluidRefValue = 0.0;
setRefCell(p, simple.dict(), pFluidRefCell, pFluidRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

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


Info<< "Reading diffusivity DTSolid\n" << endl;

dimensionedScalar DTSolid
(
    transportProperties.lookup("DTSolid")
);
dimensionedScalar DTGas
(
    transportProperties.lookup("DTGas")
);

dimensionedScalar heffGas
(
    transportProperties.lookup("heffGas")
);
dimensionedScalar heffSolid
(
    transportProperties.lookup("heffSolid")
);
//#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
//autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
//psiReactionThermo& thermo = pThermo();
//thermo.validate(args.executable(), "h", "e");


//basicSpecieMixture& composition = thermo.composition();


wordList phaseNames
    (
        transportProperties.lookup("phases")
    );

word phase1 = phaseNames[0];
word phase2 = phaseNames[1];
    
wordList speciesNames
    (
        transportProperties.lookup("speciesNames")
    ); 
    speciesTable s(wordList speciesName);

 /* wordList solidNames =speciesNames;
 wordList gasNames  =speciesNames;

 forAll(solidNames,i)
    {
        solidNames[i]+="-solid";
       
    };
    forAll(gasNames,i)
    {
         gasNames[i] +="-gas";
       
    };*/
wordList solidNames
    (
        transportProperties.lookup("solidNames")
    );
wordList gasNames
    (
        transportProperties.lookup("gasNames")
    );
    
    
basicMultiComponentMixture compositiona
    (
        transportProperties,
        solidNames,
        mesh,
        phase1
    );

basicMultiComponentMixture compositionb
    (
        transportProperties,
        gasNames,
        mesh,      
        phase2
 
    );

PtrList<volScalarField>& YSolid = compositiona.Y();
PtrList<volScalarField>& YGas = compositionb.Y();
const word inertSpecieGas(transportProperties.lookup("inertSpecieGas"));
const word inertSpecieSolid(transportProperties.lookup("inertSpecieSolid"));
if (!compositiona.species().found(inertSpecieGas))
{
    FatalIOErrorIn(args.executable().c_str(), transportProperties)
        << "Inert specie " << inertSpecieGas << " not found in available species "
        << compositiona.species() << exit(FatalIOError);
}

if (!compositionb.species().found(inertSpecieSolid))
{
    FatalIOErrorIn(args.executable().c_str(), transportProperties)
        << "Inert specie " << inertSpecieSolid << " not found in available species "
        << compositionb.species() << exit(FatalIOError);
}

dimensionedScalar massTransGas
(  
     transportProperties.lookup("massTransGas")
);



#include "createFvOptions.H"

//#include "createMRF.H"
