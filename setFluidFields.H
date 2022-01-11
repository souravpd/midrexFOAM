label pFluidRefCell = 0;
scalar pFluidRefValue = 0.0;
setRefCell(p, simple.dict(), pFluidRefCell, pFluidRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);