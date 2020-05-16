#include "chemistrySolution.H"

namespace Foam {

Ostream& operator<<(Ostream& os, const chemistrySolution& l) {
    os << l.c;
    os << l.RR;
    os << l.cellid;
    os << l.deltaTChem;
    return os;
}
Istream& operator>>(Istream& is, chemistrySolution& l) {
    is >> l.c;
    is >> l.RR;
    is >> l.cellid;
    is >> l.deltaTChem;
    return is;
}

} // namespace Foam
