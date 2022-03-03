{
    volScalarField Yt(0.0*Fe);
    // ==============Fe==============
    Info << "Solving for Fe" << endl;
    fvScalarMatrix FeSolidEqn
    (
      mvConvection->fvmDiv(phi, Fe)
      ==
      -amountConsumed
    );

    FeSolidEqn.relax();
    fvOptions.constrain(FeSolidEqn);
    FeSolidEqn.solve("Fe");

    fvOptions.correct(Fe);
    Fe.max(0.0);
    Yt = Yt + Fe;

    // =====================FeO================
    Info << "Solving for FeO" << endl;
    fvScalarMatrix FeOSolidEqn
    (
      mvConvection->fvmDiv(phi, FeO)
      ==
      amountConsumed
    );

    FeOSolidEqn.relax();
    fvOptions.constrain(FeOSolidEqn);
    FeOSolidEqn.solve("FeO");
    fvOptions.correct(FeO);
    FeO.max(0.0);
    Yt = Yt + FeO;
}
