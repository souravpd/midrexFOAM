/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "solidComponent.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidComponent::solidComponent
(
    const dictionary& solidComponentDict,
    const fvMesh& mesh
)
:
    volScalarField
    (
        IOobject
        (
            "solidC_" + solidComponentDict.dictName(),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    name_("solidC_" + solidComponentDict.dictName()),
    solidComponentDict_(solidComponentDict),
    D_(solidComponentDict_.lookup("D")),
    diffusionModel_(solidComponentDict_.lookup("diffusionModel"))
{
    Info<< "Reading field " << name_ <<"\n" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidComponent> Foam::solidComponent::clone() const
{
    return autoPtr<solidComponent>(new solidComponent(*this));
}

bool Foam::solidComponent::read(const dictionary& solidComponentDict)
{
    solidComponentDict_ = solidComponentDict;
    solidComponentDict_.lookup("diffusionModel") >> diffusionModel_;

    if (diffusionModel_ == "constant")
    {
        solidComponentDict_.lookup("D") >> D_;
        return true;
    }
    else
    {
        return false;
    }

}


// ************************************************************************* //
