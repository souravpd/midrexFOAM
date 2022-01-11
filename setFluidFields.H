label pFluidRefCell = 0;
scalar pFluidRefValue = 0.0;
setRefCell(p, simple.dict(), pFluidRefCell, pFluidRefValue);
mesh.setFluxRequired(p.name());

Info<< "Reading thermophysical properties \n" << endl;
autoPtr<psiReactionThermo> pThermoFluid(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermoFluid();
thermo.validate(args.executable() , "h" , "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if(!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
    << "Inert Specie " << inertSpecie << "not found in  available species "
    << composition.species() << exit(FatalIOError);
}

Info<< "Creating turbulence model. \n" << nl;

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> reaction
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);


Info<< "Creating field dpdtFluid\n" << endl;

volScalarField dpdtFluid
(
    IOobject
    (
        "dpdtFluid",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime , 0)
);

Info <<"Creating field kinetic energy KFluid \n" << endl;
volScalarField Kfluid("KFluid" , 0.5*magSqr(U));

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y , i)
{
    fields.add(Y[i]);
}

fields.add(thermo.he());