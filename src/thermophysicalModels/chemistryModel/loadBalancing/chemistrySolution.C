#include "chemistrySolution.H"

namespace Foam {

Ostream& operator<<(Ostream& os, const chemistrySolution& l) {
    os << l.c_increment;
    os << l.deltaTChem;
    os << l.cpuTime;
    os << l.cellid;
    return os;
}
Istream& operator>>(Istream& is, chemistrySolution& l) {
    is >> l.c_increment;
    is >> l.deltaTChem;
    is >> l.cpuTime;
    is >> l.cellid;
    return is;
}

} // namespace Foam