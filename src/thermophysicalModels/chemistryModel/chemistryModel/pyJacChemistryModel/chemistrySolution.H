#pragma once

#include "volFields.H"

namespace Foam {

struct chemistrySolution {
    scalar val;
    int    cellid;

    // TODO: implement!
    bool operator==(const chemistrySolution& rhs) const { return false; }

    bool operator!=(const chemistrySolution& rhs) const { return !(*this == rhs); }
};

// These are required for sending and receving using the OF MPI routines
Ostream& operator<<(Ostream& os, const chemistrySolution& l);
Istream& operator>>(Istream& is, chemistrySolution& l);

} // namespace Foam