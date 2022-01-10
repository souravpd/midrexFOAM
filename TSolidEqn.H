fvScalarMatrix TSolidEqn
(
   fvm::div(phiSolid, TSolid)
   - fvm::laplacian(DTSolid, TSolid) 
   + heffSolid* TSolid
   ==
    // fvOptions(TSolid) 
      heffSolid*T
);

   TSolidEqn.relax();
   fvOptions.constrain(TSolidEqn);
   TSolidEqn.solve();
   fvOptions.correct(TSolid);
