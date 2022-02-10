//Create Fluid Fields
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

dimensionedScalar DTGas
(
    transportProperties.lookup("DTGas")
);

dimensionedScalar heffGas
(
    transportProperties.lookup("heffGas")
);

//Reading Individual MassTrans Constants
dimensionedScalar massTransGasCO
(  
     transportProperties.lookup("massTransGasCO")
);
dimensionedScalar massTransGasCO2
(  
     transportProperties.lookup("massTransGasCO2")
);
dimensionedScalar massTransGasH2
(  
     transportProperties.lookup("massTransGasH2")
);
dimensionedScalar massTransGasH2O
(  
     transportProperties.lookup("massTransGasH2O")
);
dimensionedScalar massTransGasCH4
(  
     transportProperties.lookup("massTransGasCH4")
);
// //Creating massTransArray
// PtrList<dimensionedScalar> massTransGas;
// massTransGas.append(massTransGasCO);
// massTransGas.append(massTransGasCO2);
// massTransGas.append(massTransGasH2);
// massTransGas.append(massTransGasH2O);
// massTransGas.append(massTransGasCH4);
//READING INDIVIDUAL MASS FRACTIONS FOR THE SPECIES

Info << "Reading CO \n" ;
volScalarField CO
(
    IOobject 
    (
        "CO",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading CO2 \n" ;
volScalarField CO2
(
    IOobject 
    (
        "CO2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading H2 \n" ;
volScalarField H2
(
    IOobject 
    (
        "H2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading H2O \n" ;
volScalarField H2O
(
    IOobject 
    (
        "H2O",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading CH4 \n" ;
volScalarField CH4
(
    IOobject 
    (
        "CH4",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


PtrList<volScalarField> tmp_ptr_list;
tmp_ptr_list.append(CO);
tmp_ptr_list.append(CO2);
tmp_ptr_list.append(H2);
tmp_ptr_list.append(H2O);
tmp_ptr_list.append(CH4);


PtrList<volScalarField> &YGas = tmp_ptr_list;

multivariateSurfaceInterpolationScheme<scalar>:: fieldTable fields;

forAll(YGas , i){
    Info<<"Composition of the gas is" << YGas[i] << endl;
    fields.add(YGas[i]);
}

// fields.add(thermo.he());

//Create Solid Fields
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

//Set Solid Fields
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


//Set Fluid Fields

label pFluidRefCell = 0;
scalar pFluidRefValue = 0.0;
setRefCell(p, simple.dict(), pFluidRefCell, pFluidRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createFvOptions.H"