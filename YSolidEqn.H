{
  //  reaction->correct();
    volScalarField Yt(0.0*YSolid[0]);

    forAll(YSolid, i)
    {
            volScalarField& YiSolid = YSolid[i];

            fvScalarMatrix YiSolidEqn
            (
               mvConvection->fvmDiv(phi, YiSolid)
              
            // ==
              //  reaction->R(YiSolid)
              
            );

            YiSolidEqn.relax();

            fvOptions.constrain(YiSolidEqn);

            YiSolidEqn.solve("YiSolid");

            fvOptions.correct(YiSolid);

            YiSolid.max(0.0);
            Yt += YiSolid;
        
    }

    //YSolid[inertIndex] = scalar(1) - Yt;
    //YSolid[inertIndex].max(0.0);
}
