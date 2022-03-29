/*tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
   // reaction->correct();
    volScalarField Yt(0.0*YGas[0]);

    forAll(YGas, i)
    {
        if (i != inertIndexGas && compositiona.active(i))
        {
            volScalarField& YiGas = YGas[i];

            fvScalarMatrix YiGasEqn
            (
                
               mvConvection->fvmDiv(phi, YiGas)
              - fvm::laplacian(turbulence->muEff(), YiGas)
             ==
                //reaction->R(YiGas)
              fvOptions(YiGas)
            );

            YiGasEqn.relax();

            fvOptions.constrain(YiGasEqn);

            YiGasEqn.solve("YiGas");

            fvOptions.correct(YiGas);

            YiGas.max(0.0);
            Yt += YiGas;
        }
    }

    YGas[inertIndexGas] = scalar(1) - Yt;
    YGas[inertIndexGas].max(0.0);
}

*/

 dictionary gasComponentSolveDict = mesh.solutionDict().subDict("solvers").subDict("gasC_gasComponents");
        forAll(gasComponentNameList, i)
        {
            fvScalarMatrix gasMassTransport
            (
                 
                fvm::div(phi, gasComponentList[i], "div(phi,C_gasComponents)")
                - fvm::laplacian(gasComponentList[i].D(), gasComponentList[i], "laplacian(D,C_gasComponents)")
            );

            gasMassTransport.solve(gasComponentSolveDict);
        }

