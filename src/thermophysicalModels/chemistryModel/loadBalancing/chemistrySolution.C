#include "chemistrySolution.H"

namespace Foam {

Ostream& operator<<(Ostream& os, const chemistrySolution& l) {
    os << l.c_increment;
    os << l.deltaTChem;
    os << l.cpuTime;
    os << l.cellid;
    os << l.rhoi;
    return os;
}
Istream& operator>>(Istream& is, chemistrySolution& l) {
    is >> l.c_increment;
    is >> l.deltaTChem;
    is >> l.cpuTime;
    is >> l.cellid;
    is >> l.rhoi;
    return is;
}

} // namespace Foam