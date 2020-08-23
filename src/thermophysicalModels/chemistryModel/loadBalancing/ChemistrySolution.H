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
    Foam::ChemistrySolution
Description
    A small object containing everything required for updating the reaction rate
    and the chemistry time step. These are passed around in the load balancer.

SourceFiles
    ChemistrySolution.C
\*---------------------------------------------------------------------------*/

#ifndef ChemistrySolution_H
#define ChemistrySolution_H

#include "volFields.H"

namespace Foam {

struct ChemistrySolution {

    ChemistrySolution() = default;

    ChemistrySolution(label nspecie) : c_increment(nspecie, 0.0), deltaTChem(0.0), cellid(0), rhoi(0.0) {}


    // TODO: implement!
    bool operator==(const ChemistrySolution& rhs) const { return false; }

    bool operator!=(const ChemistrySolution& rhs) const { return !(*this == rhs); }

    scalarField c_increment; // =   (c_{i+1} - c_{i}) / deltaT;
    scalar deltaTChem;
    scalar cpuTime;
    label    cellid;
    scalar  rhoi;


};

///
///@brief Serialization for send
///
///@param os Buffer to serialize to
///@param s The solution to serialize
///@return Ostream& A buffer with s
///
static inline Ostream& operator<<(Ostream& os, const ChemistrySolution& s) {
    os << s.c_increment;
    os << s.deltaTChem;
    os << s.cpuTime;
    os << s.cellid;
    os << s.rhoi;
    return os;
}

///
///@brief Gets a serialized solution from Istream
///
///@param is The buffer to get the object from
///@param s  The solution to store the data
///@return Istream& A buffer without s
///
static inline Istream& operator>>(Istream& is, ChemistrySolution& s) {
    is >> s.c_increment;
    is >> s.deltaTChem;
    is >> s.cpuTime;
    is >> s.cellid;
    is >> s.rhoi;
    return is;
}


} // namespace Foam

#endif