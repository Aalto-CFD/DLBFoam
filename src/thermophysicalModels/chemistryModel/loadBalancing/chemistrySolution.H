#pragma once

#include "volFields.H"

namespace Foam {

struct chemistrySolution {

    chemistrySolution() = default;

    chemistrySolution(int nspecie) : c_increment(nspecie, 0.0), deltaTChem(0.0), cellid(0), rhoi(0.0) {}


    // TODO: implement!
    bool operator==(const chemistrySolution& rhs) const { return false; }

    bool operator!=(const chemistrySolution& rhs) const { return !(*this == rhs); }

    scalarField c_increment; // =   (c_{i+1} - c_{i}) / deltaT;
    scalar deltaTChem;
    scalar cpuTime;
    label    cellid;
    scalar  rhoi;


};

// These are required for sending and receving using the OF MPI routines
Ostream& operator<<(Ostream& os, const chemistrySolution& l);
Istream& operator>>(Istream& is, chemistrySolution& l);

} // namespace Foam