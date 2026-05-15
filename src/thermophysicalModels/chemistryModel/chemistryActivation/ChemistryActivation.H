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

Class
    Foam::ChemistryActivation

Description
    Helper routines for chemistry activation timing:
      - userTimeToTimeValue(userTime, repeat)
          Fold a user/global time into the periodic timeValue when repeat > 0.
      - chemistryActive(startTime, duration, repeat, timeValue)
          Decide whether chemistry should be active at the given timeValue,
          handling non-repeating, repeating (including wrapped intervals),
          and sentinel 'great' values indicating undefined options.
\*---------------------------------------------------------------------------*/

#ifndef ChemistryActivation_H
#define ChemistryActivation_H

#include "OSspecific.H"

namespace Foam
{

inline bool chemistryActive
(
    const scalar startTime,
    const scalar duration,
    const scalar repeat,
    const scalar timeValue
)
{
    if (startTime <= -great) return true;
    if (duration <= scalar(0)) return false;

    const scalar time = timeValue;  // keep global time

    // non-repeating case
    if (repeat <= 0)
    {
        return (time >= startTime && time <= startTime + duration);
    }

    // before first activation
    if (time < startTime) return false;

    // infinite active window
    if (duration >= repeat) return true;

    const scalar dt = time - startTime;
    const scalar modt = fmod(dt, repeat);

    return (modt >= 0 && modt <= duration);
}

} // namespace Foam

#endif