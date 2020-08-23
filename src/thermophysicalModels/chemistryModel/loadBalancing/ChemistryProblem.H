/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::ChemistryProblem
Description
    A small object containing everything required for solving the reaction rate
    using the ODE. These are passed around in the load balancer.

SourceFiles
    ChemistryProblem.C
\*---------------------------------------------------------------------------*/

#ifndef ChemistryProblem_H
#define ChemistryProblem_H

#include "volFields.H"

namespace Foam
{

struct ChemistryProblem
{

    ChemistryProblem() = default;
    ChemistryProblem(label nSpecie)
        : c(nSpecie), Ti(0), pi(0), rhoi(0), deltaTChem(0), deltaT(0), cellid(0)
    {
    }

    scalarField c;
    scalar      Ti;
    scalar      pi;
    scalar      rhoi;
    scalar      deltaTChem;
    scalar      deltaT;
    scalar      cpuTime;
    label       cellid;

    // TODO: implement!
    bool operator==(const ChemistryProblem& rhs) const
    {
        return false;
    }

    bool operator!=(const ChemistryProblem& rhs) const
    {
        return !(*this == rhs);
    }
};

///
///@brief Serialization for send
///
///@param os Buffer to serialize to
///@param p Problem to serialize
///@return Ostream& A buffer with p
///
static inline Ostream& operator<<(Ostream& os, const ChemistryProblem& p)
{

    os << p.c;
    os << p.Ti;
    os << p.pi;
    os << p.rhoi;
    os << p.deltaTChem;
    os << p.deltaT;
    os << p.cpuTime;
    os << p.cellid;

    return os;
}

///
///@brief Gets a serialized problem from Istream
///
///@param is The buffer to get the object from
///@param p  The problem to store the data
///@return Istream& A buffer without p
///
static inline Istream& operator>>(Istream& is, ChemistryProblem& p)
{

    is >> p.c;
    is >> p.Ti;
    is >> p.pi;
    is >> p.rhoi;
    is >> p.deltaTChem;
    is >> p.deltaT;
    is >> p.cpuTime;
    is >> p.cellid;

    return is;
}

} // namespace Foam

#endif