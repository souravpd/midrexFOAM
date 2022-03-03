/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    midrexFoam

Description
    Midrex flow solver which solves for the solid velocity potential, to
    calculate the flux-field, from which the velocity field is obtained by
    reconstructing the flux.

    This application is particularly useful to generate starting fields for
    Navier-Stokes codes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "nonOrthogonalSolutionControl.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "multiComponentMixture.H"
#include "psiReactionThermo.H"
#include "psiThermo.H"
#include "CombustionModel.H"
#include "turbulentFluidThermoModel.H"
#include "fixedGradientFvPatchFields.H"
#include "pressureControl.H"
#include "coordinateSystem.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "multivariateScheme.H"
#include "fvcSmooth.H"
#include "localEulerDdtScheme.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "pName",
        "pName",
        "Name of the pressure field"
    );

    argList::addBoolOption
    (
        "initialiseUBCs",
        "Initialise U boundary conditions"
    );

    argList::addBoolOption
    (
        "writePhi",
        "Write the velocity potential field"
    );

    argList::addBoolOption
    (
        "writep",
        "Calculate and write the pressure field"
    );

    argList::addBoolOption
    (
        "withFunctionObjects",
        "execute functionObjects"
    );

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "initContinuityErrs.H"


    simpleControl simple(mesh, "simple");

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< nl << "Calculating potential flow" << endl;

    // Since solver contains no time loop it would never execute
    // function objects so do it ourselves
    runTime.functionObjects().start();

    // Non-orthogonal velocity potential corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix PhiSolidEqn
        (
            fvm::laplacian(dimensionedScalar(dimless, 1), PhiSolid)
         ==
            fvc::div(phiSolid)
        );

        PhiSolidEqn.setReference(PhiRefCell, PhiRefValue);
        PhiSolidEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phiSolid -= PhiSolidEqn.flux();
        }
    }


    //MRF.makeAbsolute(phi);

    Info<< "Continuity error = "
        << mag(fvc::div(phiSolid))().weightedAverage(mesh.V()).value()
        << endl;

    USolid = fvc::reconstruct(phiSolid);
    USolid.correctBoundaryConditions();
 

      Info<< "Interpolated velocity error = "
        << (sqrt(sum(sqr(fvc::flux(USolid) - phiSolid)))/sum(mesh.magSf())).value()
        << endl;
   
    Info<< nl << "\nWriting U and phi for solid\n" << endl;
    USolid.write();
    phiSolid.write();
    
    // Optionally write Phi
    if (args.optionFound("writePhi"))
    {
        PhiSolid.write();
    }
    runTime.functionObjects().end();
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    p.storePrevIter();

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
            #include "TEqn.H"
            #include "TSolidEqn.H"
            #include "YGasEqn.H"
            #include "YSolidEqn.H"
        }

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();
       
        T.write();
        TSolid.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
