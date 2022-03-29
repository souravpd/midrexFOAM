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
    Foam::solidComponent

Description
    Single incompressible solution solidComponent.
    Used as part of the multiComponentSolution to store solidComponent data.

SourceFiles
    solidComponent.C

Author
    Rudolf Hellmuth (rudolf.hellmuth@gmail.com)

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "dictionaryEntry.H"

#ifndef solidComponent_H
#define solidComponent_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidComponent Declaration
\*---------------------------------------------------------------------------*/

class solidComponent
:
    public volScalarField
{
    // Private data

        word name_;
        dictionary solidComponentDict_;

        dimensionedScalar D_;
        word diffusionModel_;

public:

    // Constructors

        //- Construct from components
        solidComponent
        (
            const dictionary& solidComponentDict,
            const fvMesh& mesh
        );

        //- Return clone
        autoPtr<solidComponent> clone() const;

        //- Return a pointer to a new solidComponent created on freestore
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
            autoPtr<solidComponent> operator()(Istream& is) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return autoPtr<solidComponent>(new solidComponent(ent, mesh_));
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

        //- Return const-access to solidComponent diffusion coefficient
        const dimensionedScalar& D() const
        {
            return D_;
        }

        //- Return const-access to solidComponent diffusion model
        const word& diffusionModel() const
        {
            return diffusionModel_;
        }

        //- Return const-access to solidComponent molar mass
//        const dimensionedScalar& molarMass() const
//        {
//            return molarMass_;
//        }

        //- Read base transportProperties dictionary
        bool read(const dictionary& solidComponentDict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
