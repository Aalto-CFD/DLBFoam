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
#line 59 "/u/97/tekgulb1/data/Documents/FOAM-Aalto/tutorials/reactingFoam/shearlayer/0/U/#codeStream"
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
    void codeStream_857e5d87771f7ebb2d48672cf7fdabfb1060b378
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 24 "/u/97/tekgulb1/data/Documents/FOAM-Aalto/tutorials/reactingFoam/shearlayer/0/U/#codeStream"
const IOdictionary& d = static_cast<const IOdictionary&>(dict);
        const fvMesh& mesh = refCast<const fvMesh>(d.db());
        const vectorField& CC = mesh.C(); //cell center
        
         
        vectorField U(mesh.nCells());
        scalar Umax = 40;

        scalar A = 500;
        scalar L = 0.008;
        scalar ro = 0.0004;
        scalar B = 5;
        
        

        forAll(CC,cellI)
        {

            
        	scalar x = CC[cellI].x();
        	scalar y = CC[cellI].y(); 
        	scalar z = CC[cellI].z(); 	

            U[cellI] = vector(1,0,0) * Umax * 0.5 *(1-tanh(B*(mag(y-L/2)/ro - ro/mag(y-L/2)))) + vector(0,1,0)*(Umax/2000)*sin(2*constant::mathematical::pi*x/(L/2));
            
    
        }
        
        writeEntry(os,"",U);
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

