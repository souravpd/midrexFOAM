tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

dimensionedScalar massTransGas[5];
massTransGas[0] = massTransGasCO;
massTransGas[1] = massTransGasCO2;
massTransGas[2] = massTransGasH2;
massTransGas[3] = massTransGasH2O;
massTransGas[4] = massTransGasCH4;
{
   // reaction->correct();
    volScalarField Yt(0.0*YGas[0]);

    forAll(YGas, i)
    {


            Info << "solving YGas" << YGas[i] << nl;
            Info << "massTransGas" << massTransGas[i] << nl;
            volScalarField& YiGas = YGas[i];

            fvScalarMatrix YiGasEqn
            (
              mvConvection->fvmDiv(phi, YiGas)
              ==
              fvm::laplacian(massTransGas[i], YiGas) 
            );

            YiGasEqn.relax();

            fvOptions.constrain(YiGasEqn);

            YiGasEqn.solve("YiGas");

            fvOptions.correct(YiGas);

            YiGas.max(0.0);
            Yt += YiGas;
        
    }

   // YGas[inertIndex] = scalar(1) - Yt;
    // YGas[inertIndex].max(0.0);
}
