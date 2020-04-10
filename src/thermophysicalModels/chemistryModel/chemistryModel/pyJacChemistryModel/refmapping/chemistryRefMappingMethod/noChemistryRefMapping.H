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
    Foam::chemistryRefMappingMethods::none

Description

\*---------------------------------------------------------------------------*/

#ifndef noChemistryRefMapping_H
#define noChemistryRefMapping_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chemistryRefMappingMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class none Declaration
\*---------------------------------------------------------------------------*/


///
///@brief Not sure if this class is necessary because theres the switch "active" on/off in the base class.
///       Maybe this can be used for debugging at some point
///
///@tparam CompType NO IDEA
///@tparam ThermoType NO IDEA
///
template<class CompType, class ThermoType>
class noChemistryRefMapping : public chemistryRefMappingMethod<CompType, ThermoType>{

public:

    //Note this is defined here empty.
    virtual void applyMapping(pyJacChemistryModel<CompType, ThermoType>& chemistry) const {};

protected:

    // None

private:

    // None


};


    


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//No idea what this means so always include
//#ifdef NoRepository
    #include "noChemistryRefMapping.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
