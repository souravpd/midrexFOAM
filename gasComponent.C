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

#include "gasComponent.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gasComponent::gasComponent
(
    const dictionary& gasComponentDict,
    const fvMesh& mesh
)
:
    volScalarField
    (
        IOobject
        (
            "gasC_" + gasComponentDict.dictName(),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    name_("gasC_" + gasComponentDict.dictName()),
    gasComponentDict_(gasComponentDict),
    D_(gasComponentDict_.lookup("D")),
    diffusionModel_(gasComponentDict_.lookup("diffusionModel"))
{
    Info<< "Reading field " << name_ <<"\n" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::gasComponent> Foam::gasComponent::clone() const
{
    return autoPtr<gasComponent>(new gasComponent(*this));
}

bool Foam::gasComponent::read(const dictionary& gasComponentDict)
{
    gasComponentDict_ = gasComponentDict;
    gasComponentDict_.lookup("diffusionModel") >> diffusionModel_;

    if (diffusionModel_ == "constant")
    {
        gasComponentDict_.lookup("D") >> D_;
        return true;
    }
    else
    {
        return false;
    }

}


// ************************************************************************* //
