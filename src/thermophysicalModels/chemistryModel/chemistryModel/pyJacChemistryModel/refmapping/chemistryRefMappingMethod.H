/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::chemistryRefMappingMethod

Description
    An abstract class for chemistry tabulation.

SourceFiles
    chemistryRefMappingMethod.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryRefMappingMethod_H
#define chemistryRefMappingMethod_H

#include "IOdictionary.H"
#include "scalarField.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"



// Mixture fraction headers
#include "atomicWeights.H"
#include "psiReactionThermo.H"
#include "mixtureFraction.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                           Class chemistryRefMappingMethod Declaration
\*---------------------------------------------------------------------------*/

class chemistryRefMappingMethod
{

public:


    // Constructors
    chemistryRefMappingMethod
    (
        const dictionary& dict,
        const wordList& species
    );

    ///
    ///@brief THIS IS THE BREAD AND BUTTER WHICH SHOULD BE DEFINED IN ALL CLASSES!
    ///       THIS IS NOT PURE AT THIS POINT 
    ///
    // TODO add params
    virtual void applyMapping() const = 0;


    ///
    ///@brief Dont know.
    ///
    ///@param composition 
    ///
    //this is a bit useless and should be moved to constructor
    void init_mixture_fraction(basicSpecieMixture& composition);
    

    ///
    ///@brief Checks if a reference cell or not
    ///
    ///@param Y 
    ///@param celli 
    ///@return true 
    ///@return false 
    ///
    bool check_if_refcell(PtrList<volScalarField>& Y, const label celli);

    
    


    //- Destructor
    virtual ~chemistryRefMappingMethod() {}


    
    ///
    ///@brief Is the reference cell mapping active or not
    ///
    ///@return true 
    ///@return false 
    ///
    inline bool active()
    {
        return active_;
    }
    


protected:

    const dictionary& dict_;

    const dictionary coeffsDict_;

    //- Is tabulation active?
    Switch active_;

    int mytest;
    mixtureFraction mixture_fraction_;


private:


};



} //namespace Foam


#endif

// ************************************************************************* //
