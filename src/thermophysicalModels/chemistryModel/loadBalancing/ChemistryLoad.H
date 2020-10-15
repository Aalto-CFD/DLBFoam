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
    Foam::ChemistryLoad

Description
    A generic load object that stores a rank and a value and enables comparison
as well as serialization through the ostream and istream operators.

SourceFiles
    ChemistryLoad.C
\*---------------------------------------------------------------------------*/

#ifndef ChemistryLoad_H
#define ChemistryLoad_H

#include "Istream.H"
#include "Ostream.H"

namespace Foam
{

struct ChemistryLoad
{
    label  rank;  // rank of the process
    scalar value; // an indicator how busy the processor is

    ChemistryLoad() = default;

    ChemistryLoad(label _rank, scalar _value) : rank(_rank), value(_value)
    {
    }

    // overload comparisons for sorting
    bool operator==(const ChemistryLoad& rhs) const
    {
        return value == rhs.value;
    }
    bool operator!=(const ChemistryLoad& rhs) const
    {
        return !(value == rhs.value);
    }
    bool operator<=(const ChemistryLoad& rhs) const
    {
        return value <= rhs.value;
    }
    bool operator>=(const ChemistryLoad& rhs) const
    {
        return value >= rhs.value;
    }
    bool operator<(const ChemistryLoad& rhs) const
    {
        return value < rhs.value;
    }
    bool operator>(const ChemistryLoad& rhs) const
    {
        return value > rhs.value;
    }
};
///
///@brief Serialization for send ostream
///
///@param os Buffer to serialize to
///@param l  Load to serialize
///@return Ostream& Buffer with the object l
///
static inline Ostream& operator<<(Ostream& os, const ChemistryLoad& l)
{

    os << l.rank;
    os << l.value;

    return os;
}

///
///@brief Gets a serialized load from istream
///
///@param is Buffer to get from
///@param l  Load to which to assign
///@return Istream& Buffer without l
///
static inline Istream& operator>>(Istream& is, ChemistryLoad& l)
{

    is >> l.rank;
    is >> l.value;

    return is;
}

} // namespace Foam

#endif