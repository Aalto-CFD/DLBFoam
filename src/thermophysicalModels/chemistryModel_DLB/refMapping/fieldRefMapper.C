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

#include "fieldRefMapper.H"

Foam::fieldRefMapper::fieldRefMapper(const dictionary& dict)
:
    dict_(dict),
    active_(dict_.lookupOrDefault<Switch>("active", false)),
    fieldName_(word::null),
    minValue_(-VGREAT),
    maxValue_(VGREAT),
    Ttolerance_(dict_.lookupOrDefault<scalar>("deltaT", VGREAT))
{
    if (active_)
    {
        fieldName_ = dict_.lookup<word>("field");
        minValue_ = dict_.lookup<scalar>("min");
        maxValue_ = dict_.lookup<scalar>("max");

        if (minValue_ > maxValue_)
        {
            FatalIOErrorInFunction(dict_)
                << "Minimum field value " << minValue_
                << " exceeds maximum field value " << maxValue_
                << exit(FatalIOError);
        }
    }
}


bool Foam::fieldRefMapper::shouldMap(const scalar fieldValue) const
{
    return active_ && fieldValue >= minValue_ && fieldValue <= maxValue_;
}


bool Foam::fieldRefMapper::temperatureWithinRange
(
    const scalar Ti,
    const scalar Tref
) const
{
    return abs(Ti - Tref) < Ttolerance_;
}


// ************************************************************************* //