/*---------------------------------------------------------------------------*  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 60 "/u/97/tekgulb1/data/Documents/FOAM-Aalto/tutorials/reactingFoam/shearlayer/0/NC12H26/#codeStream"
#include "fvCFD.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    void codeStream_94c926906ae2f61c6744a8173c211f9cb772ae65
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 25 "/u/97/tekgulb1/data/Documents/FOAM-Aalto/tutorials/reactingFoam/shearlayer/0/NC12H26/#codeStream"
const IOdictionary& d = static_cast<const IOdictionary&>(dict);
        const fvMesh& mesh = refCast<const fvMesh>(d.db());
        const vectorField& CC = mesh.C(); //cell center
        
         
        scalarField NC12H26(mesh.nCells());
        scalar NC12H26_i = 0.0;

        scalar A = 1;
        scalar L = 0.008;
        scalar ro = 0.0004;
        scalar B = 5;
        
        

        forAll(CC,cellI)
        {

            
        	scalar x = CC[cellI].x();
        	scalar y = CC[cellI].y(); 
        	scalar z = CC[cellI].z(); 	

            
            NC12H26[cellI] = NC12H26_i + A * 0.5 * (1-tanh(B*(mag(y-L/2)/ro - ro/mag(y-L/2))));
            
        }
        
        writeEntry(os,"",NC12H26);
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

