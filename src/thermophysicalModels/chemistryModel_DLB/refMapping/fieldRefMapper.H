/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           |
     \\/     M anipulation  | 2020, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of DLBFoam library, derived from OpenFOAM.

    https://github.com/Aalto-CFD/DLBFoam

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

#ifndef fieldRefMapper_H
#define fieldRefMapper_H

#include "IOdictionary.H"
#include "Switch.H"

namespace Foam
{

class fieldRefMapper
{
public:

    fieldRefMapper() = default;

    fieldRefMapper(const dictionary& dict);

    //- Check if the field value is within the configured range
    bool shouldMap(const scalar fieldValue) const;

    //- Check if the temperature is within the configured tolerance
    bool temperatureWithinRange(const scalar Ti, const scalar Tref) const;

    //- Is reference mapping active?
    bool active() const
    {
        return active_;
    }

    //- Name of the field used to select cells for mapping
    const word& fieldName() const
    {
        return fieldName_;
    }


private:

    const dictionary dict_;

    Switch active_;

    word fieldName_;

    scalar minValue_;

    scalar maxValue_;

    //- Tolerance for temperature checking [K]
    scalar Ttolerance_;
};

} // namespace Foam

#endif

// ************************************************************************* //