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