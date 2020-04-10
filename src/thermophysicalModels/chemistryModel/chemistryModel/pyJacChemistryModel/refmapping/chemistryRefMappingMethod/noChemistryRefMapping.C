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

\*---------------------------------------------------------------------------*/

#include "noChemistryRefMapping.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
namespace chemistryRefMappingMethods
{


template<class CompType, class ThermoType>
noChemistryRefMapping<CompType, ThermoType>::noChemistryRefMapping
(
    const dictionary& chemistryProperties,
    pyJacChemistryModel<CompType, ThermoType>& chemistry
)
:
    chemistryRefMappingMethod<CompType, ThermoType>
    (
        chemistryProperties,
        chemistry
    )
{
    this->active_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
noChemistryRefMapping<CompType, ThermoType>::~noChemistryRefMapping()
{}

//Will be optimized away...
template<class CompType, class ThermoType>
noChemistryRefMapping<CompType, ThermoType>::~applyMapping()
{}


} //namespace Foam

} //namespace chemistryRefMappingMethods

// ************************************************************************* //
