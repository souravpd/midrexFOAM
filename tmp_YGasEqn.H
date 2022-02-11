tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);
{
   // reaction->correct();
    volScalarField Yt(0.0*YGas[0]);

    forAll(YGas, i)
    {


            // Info << "solving YGas" <<  nl;
            volScalarField& YiGas = YGas[i];
            Info << "Ygas[i]" << YiGas << endl;
            fvScalarMatrix YiGasEqn
            (
              mvConvection->fvmDiv(phi, YiGas)
              ==
              fvm::laplacian(dimensionedScalar("massTransGas", dimensionSet(0, 2, -1, 0, 0 , 0 ,0), vals[i]), YiGas) 
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
