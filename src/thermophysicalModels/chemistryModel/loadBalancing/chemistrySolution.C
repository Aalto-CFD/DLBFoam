#include "chemistrySolution.H"

namespace Foam {

Ostream& operator<<(Ostream& os, const chemistrySolution& l) {
    os << l.c_increment;
    os << l.cellid;
    os << l.deltaTChem;
    return os;
}
Istream& operator>>(Istream& is, chemistrySolution& l) {
    is >> l.c_increment;
    is >> l.cellid;
    is >> l.deltaTChem;
    return is;
}

} // namespace Foam