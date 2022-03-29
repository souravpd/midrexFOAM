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

Class
    Foam::gasComponent

Description
    Single incompressible solution gasComponent.
    Used as part of the multiComponentSolution to store gasComponent data.

SourceFiles
    gasComponent.C

Author
    Rudolf Hellmuth (rudolf.hellmuth@gmail.com)

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "dictionaryEntry.H"

#ifndef gasComponent_H
#define gasComponent_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class gasComponent Declaration
\*---------------------------------------------------------------------------*/

class gasComponent
:
    public volScalarField
{
    // Private data

        word name_;
        dictionary gasComponentDict_;

        dimensionedScalar D_;
        word diffusionModel_;

public:

    // Constructors

        //- Construct from components
        gasComponent
        (
            const dictionary& gasComponentDict,
            const fvMesh& mesh
        );

        //- Return clone
        autoPtr<gasComponent> clone() const;

        //- Return a pointer to a new gasComponent created on freestore
        //  from Istream
        class iNew
        {
            const fvMesh& mesh_;

        public:

            iNew
            (
                const fvMesh& mesh
            )
            :
                mesh_(mesh)
            {}
            autoPtr<gasComponent> operator()(Istream& is) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return autoPtr<gasComponent>(new gasComponent(ent, mesh_));
            }
        };

    // Destructor


    // Member Functions

        const word& name() const
        {
            return name_;
        }

        const word& keyword() const
        {
            return name();
        }

        //- Return const-access to gasComponent diffusion coefficient
        const dimensionedScalar& D() const
        {
            return D_;
        }

        //- Return const-access to gasComponent diffusion model
        const word& diffusionModel() const
        {
            return diffusionModel_;
        }

        //- Return const-access to gasComponent molar mass
//        const dimensionedScalar& molarMass() const
//        {
//            return molarMass_;
//        }

        //- Read base transportProperties dictionary
        bool read(const dictionary& gasComponentDict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
