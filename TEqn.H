fvScalarMatrix TEqn
(
   fvm::div(phi, T)
   - fvm::laplacian(DTGas, T) 
   + heffGas* T
   ==
    // fvOptions(TSolid) 
       heffGas * TSolid
);

   TEqn.relax();
   fvOptions.constrain(TEqn);
   TEqn.solve();
   fvOptions.correct(T);
