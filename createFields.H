//Create Fluid Fields
#include "createRDeltaT.H"
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



Info<< "end of reading Transport Properties\n" << endl;
Info<< "Reading solid Diffusion Properties\n" << endl;
IOdictionary solidTransportProperties
(
    IOobject
    (
        "solidTransportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

 

//IOdictionary::toc() : table of contents
const wordList solidComponentNameList(solidTransportProperties.toc());

PtrList<solidComponent> solidComponentList(solidComponentNameList.size());

//for divSchemes later on
//multivariateSurfaceInterpolationScheme<scalar>::fieldTable gasComponentFields;

forAll(solidComponentNameList, i)
{
    const word& solidComponentName = solidComponentNameList[i];
    const dictionary& solidComponentDict = solidTransportProperties.subDict(solidComponentName);

    solidComponentList.set
    (
        i,
        new solidComponent
        (
            solidComponentDict,
            mesh
        )
    );


}






Info<< "Reading gas Diffusion Properties\n" << endl;
IOdictionary gasTransportProperties
(
    IOobject
    (
        "gasTransportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

 

//IOdictionary::toc() : table of contents
const wordList gasComponentNameList(gasTransportProperties.toc());

PtrList<gasComponent> gasComponentList(gasComponentNameList.size());

//for divSchemes later on
//multivariateSurfaceInterpolationScheme<scalar>::fieldTable gasComponentFields;

forAll(gasComponentNameList, i)
{
    const word& gasComponentName = gasComponentNameList[i];
    const dictionary& gasComponentDict = gasTransportProperties.subDict(gasComponentName);

    gasComponentList.set
    (
        i,
        new gasComponent
        (
            gasComponentDict,
            mesh
        )
    );
}





//    gasComponentFields.add(gasComponentList[i]);



//    Info<< "\nDiffusion = " << gasComponentList[0].name() << endl;
/*
Info<< "Reading thermophysical properties of gas\n" << endl;
autoPtr<rhoReactionThermo> pThermo(rhoReactionThermo::New(Mesh));
rhoReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& compositionGas = thermo.composition();
PtrList<volScalarField>& YGas = compositionGas.Y();

const word inertSpecieGas(thermo.lookup("inertSpecieGas"));
if (!compositiona.species().found(inertSpecieGas))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecieGas << " not found in available species "
        << compositiona.species()
        << exit(FatalIOError);
}

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(YGas, i)
{
    fields.add(YGas[i]);
}
fields.add(thermo.he());*/

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


dimensionedScalar DTGas
(
    transportProperties.lookup("DTGas")
);

dimensionedScalar heffGas
(
    transportProperties.lookup("heffGas")
);
/*
Info<< "Reading thermophysical properties of Solid\n" << endl;
basicSpecieMixture& compositionSolid = thermo.compositionb();
PtrList<volScalarField>& YSolid = compositionSolid.Y();

const word inertSpecieSolid(thermo.lookup("inertSpecieSolid"));
if (!compositionb.species().found(inertSpecieSolid))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecieSolid << " not found in available species "
        << compositionb.species()
        << exit(FatalIOError);
}


forAll(YSolid, i)
{
    fields.add(YSolid[i]);
}
fields.add(thermo.he());

*/
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
Info<< "Creating reaction model\n" << endl;
/*autoPtr<CombustionModel<rhoReactionThermo>> reaction
(
    CombustionModel<rhoReactionThermo>::New(thermo, turbulence())
);*/
#include "createFvOptions.H"



