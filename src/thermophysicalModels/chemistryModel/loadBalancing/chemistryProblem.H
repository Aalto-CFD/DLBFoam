#pragma once

#include "volFields.H"

namespace Foam {

struct chemistryProblem {

    chemistryProblem() = default;
    chemistryProblem(int nSpecie) : c(nSpecie), Ti(0), pi(0), rhoi(0), deltaTChem(0), deltaT(0), cellid(0){}

    scalarField c;
    scalar      Ti;
    scalar      pi;
    scalar      rhoi;
    scalar      deltaTChem;
    scalar      deltaT;
    scalar      cpuTime;
    label       cellid;

    // TODO: implement!
    bool operator==(const chemistryProblem& rhs) const { return false; }

    bool operator!=(const chemistryProblem& rhs) const { return !(*this == rhs); }
};

// These are required for sending and receving using the OF MPI routines
Ostream& operator<<(Ostream& os, const chemistryProblem& l);
Istream& operator>>(Istream& is, chemistryProblem& l);

} // namespace Foam