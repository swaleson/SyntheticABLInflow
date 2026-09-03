/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::synInflowFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    This boundary condition describes a synthetic inflow condition corrected by applying the VBIC procedure

Usage
    No arguments to be provided apart from the VBICDict in the system folder

Note

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types

SourceFiles
    synInflowFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef synInflowFvPatchVectorField_H
#define synInflowFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
         Class synInflowFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class synInflowFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
        #include "initializeInflowConditions.H"

        // Initializing local variables
        vector      sf_xDir;
        vector      sf_yDir;

        scalar      sf_uSmall;
        scalar      sf_lSmall;

        bool sf_calculateCorrections;
        vectorField sf_faceC;
        scalarField sf_nodalArea;
        vectorField sf_nodesTrans;
        vectorField sf_nodesProps;
        vectorField sf_seedsWeights;
        scalarField sf_faceAreas;
        scalarField sf_facesIds;
        vectorField sf_seedsProps;
        vectorField sf_uTot_c;
        scalar      sf_dxScale;
        int         sf_curTimeIndex;

        // For LU solver
        scalarField ls_Ld;
        vectorField ls_Lo;
        scalarField ls_Ud;
        vectorField ls_Uo;
        scalarField ls_LUpc ;
        scalarField ls_LUpr;

public:

   //- Runtime type information
   TypeName("synInflow");


   // Constructors

        //- Construct from patch and internal field
        synInflowFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        synInflowFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping
        synInflowFvPatchVectorField
        (
            const synInflowFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        synInflowFvPatchVectorField
        (
            const synInflowFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new synInflowFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        synInflowFvPatchVectorField
        (
            const synInflowFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new synInflowFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// ************************************************************************* //

