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

    //- Construct from components
    chemistryRefMappingMethod
    (
        const dictionary& dict,
        const wordList& species
    )
    :
    dict_(dict),
    coeffsDict_(dict.subDict("refmapping")),
    active_(coeffsDict_.lookupOrDefault<Switch>("active", false)),
    mixture_fraction_(coeffsDict_.subDict("mixtureFractionProperties"), species)
    {


        if (this->active_)
        {
        // scalar test = readScalar(mixFracDict.lookup("test"));
        }
    }


    ///
    ///@brief THIS IS THE BREAD AND BUTTER WHICH SHOULD BE DEFINED IN ALL CLASSES!
    ///       THIS IS NOT PURE AT THIS POINT 
    ///
    // TODO add params
    virtual void applyMapping() const = 0;


    

    //this is a bit useless and should be moved to constructor
    void init_mixture_fraction(basicSpecieMixture& composition){
      
        mixture_fraction_.update(composition);

    };


    bool check_if_refcell(PtrList<volScalarField>& Y, const label celli)
    {

        auto beta_of = mixture_fraction_.get_beta();
        auto alpha = mixture_fraction_.get_alpha();

        scalar beta = 0.0; //TODO: rename!
        scalar Z;
        forAll(Y, iField)
        {
            const scalarField& Yi = Y[iField];
            beta += alpha[iField]*Yi[celli];
        }
        Z = (beta - beta_of[0])/(beta_of[1] - beta_of[0]);
        if (Z>1e-5)
        {
            return false;
        }
        else
        {
            return true;
        }

    }


    //- Destructor
    virtual ~chemistryRefMappingMethod() {}


    // Member Functions

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
