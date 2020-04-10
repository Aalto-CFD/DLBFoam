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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
class pyJacChemistryModel;

/*---------------------------------------------------------------------------*\
                           Class chemistryRefMappingMethod Declaration
\*---------------------------------------------------------------------------*/

///
///@brief A base class for chemistry reference cell mappings
///
///@tparam CompType NO IDEA
///@tparam ThermoType NO IDEA
///
template<class CompType, class ThermoType>
class chemistryRefMappingMethod
{

public:

    //- Runtime type information
    //TypeName("chemistryRefMappingMethod");


    

    // Constructors

    //- Construct from components
    chemistryRefMappingMethod
    (
        const dictionary& dict
    );


    ///
    ///@brief THIS IS THE BREAD AND BUTTER WHICH SHOULD BE DEFINED IN ALL CLASSES!
    ///
    // TODO add params
    virtual void applyMapping(pyJacChemistryModel<CompType, ThermoType>& chemistry) const = 0;




    //- Destructor
    virtual ~chemistryRefMappingMethod();


    // Member Functions

    inline bool active()
    {
        return active_;
    }

    inline bool log()
    {
        return active_ && log_;
    }

    inline bool variableTimeStep()
    {
        return chemistry_.variableTimeStep();
    }

    inline scalar tolerance() const
    {
        return tolerance_;
    }

        


protected:

    const dictionary& dict_;

    const dictionary coeffsDict_;

    //- Is tabulation active?
    Switch active_;

    //- Switch to select performance logging
    Switch log_;

    //TODO I dont think this should have a chemistryModel member variable,
    //certainly not a specific pyJacModel
    pyJacChemistryModel<CompType, ThermoType>& chemistry_;

    scalar tolerance_;



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



} //namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//TODO: find out when NoRepositiory is turned off
//#ifdef NoRepository
    #include "chemistryRefMappingMethod.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
