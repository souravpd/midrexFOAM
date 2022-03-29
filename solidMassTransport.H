/*tmp<fv::convectionScheme<scalar>> mvConvectionSolid
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phiSolid,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
   // reaction->correct();
    volScalarField Yt(0.0*YSolid[0]);

    forAll(YSolid, i)
    {
        if (i != inertIndexSolid && compositionb.active(i))
        {
            volScalarField& YiSolid = YSolid[i];

            fvScalarMatrix YiSolidEqn
            (
                
               mvConvectionSolid->fvmDiv(phi, YiSolid)
             
             ==
                //reaction->R(Yisolid)
              fvOptions(YiSolid)
            );

            YiSolidEqn.relax();

            fvOptions.constrain(YiSolidEqn);

            YiSolidEqn.solve("YiSolid");

            fvOptions.correct(YiSolid);

            YiSolid.max(0.0);
            Yt += YiSolid;
        }
    }

    YSolid[inertIndexSolid] = scalar(1) - Yt;
    YSolid[inertIndexSolid].max(0.0);
}*/
dictionary solidComponentSolveDict = mesh.solutionDict().subDict("solvers").subDict("solidC_solidComponents");
        forAll(solidComponentNameList, i)
        {
            fvScalarMatrix solidMassTransport
            (
                
                fvm::div(phiSolid, solidComponentList[i], "div(phi,C_solidComponents)")
                - fvm::laplacian(solidComponentList[i].D(), solidComponentList[i], "laplacian(D,C_solidComponents)")
            );

            solidMassTransport.solve(solidComponentSolveDict);
        }

