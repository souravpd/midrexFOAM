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
    volScalarField Yt(0.0*CO);

//========================CO========================

    Info << "Solving for CO" << endl;
    fvScalarMatrix COGasEqn
    (
        mvConvection->fvmDiv(phi , CO)
        ==
        fvm::laplacian(massTransGasCO , CO)
    );

    COGasEqn.relax();
    fvOptions.constrain(COGasEqn);

    COGasEqn.solve("CO");

    fvOptions.correct(CO);
    CO.max(0.0);
    Yt = Yt + CO;

//===========================CO2=======================
    Info << "Solving for CO2" << endl;
    fvScalarMatrix CO2GasEqn
    (
        mvConvection->fvmDiv(phi , CO2)
        ==
        fvm::laplacian(massTransGasCO2 , CO2)
    );

    CO2GasEqn.relax();
    fvOptions.constrain(CO2GasEqn);

    CO2GasEqn.solve("CO2");

    fvOptions.correct(CO2);
    CO2.max(0.0);
    Yt = Yt + CO2;

//===========================H2=======================
    Info << "Solving for H2" << endl;
    fvScalarMatrix H2GasEqn
    (
        mvConvection->fvmDiv(phi , H2)
        ==
        fvm::laplacian(massTransGasH2 , H2)
    );

    H2GasEqn.relax();
    fvOptions.constrain(H2GasEqn);

    H2GasEqn.solve("H2");

    fvOptions.correct(H2);
    H2.max(0.0);
    Yt = Yt + H2;

//===========================H2O=======================
    Info << "Solving for H2O" << endl;
    fvScalarMatrix H2OGasEqn
    (
        mvConvection->fvmDiv(phi , H2O)
        ==
        fvm::laplacian(massTransGasH2O , H2O)
    );

    H2OGasEqn.relax();
    fvOptions.constrain(H2OGasEqn);

    H2OGasEqn.solve("H2O");

    fvOptions.correct(H2O);
    H2O.max(0.0);
    Yt = Yt + H2O;


//===========================CH4=======================
    Info << "Solving for CH4" << endl;
    fvScalarMatrix CH4GasEqn
    (
        mvConvection->fvmDiv(phi , CH4)
        ==
        fvm::laplacian(massTransGasCH4 , CH4)
    );

    CH4GasEqn.relax();
    fvOptions.constrain(CH4GasEqn);

    CH4GasEqn.solve("CH4");

    fvOptions.correct(CH4);
    CH4.max(0.0);
    Yt = Yt + CH4;
}

//THis code belongs to a new branch