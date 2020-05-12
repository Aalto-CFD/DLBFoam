#include "chemistrySolution.H"

namespace Foam {

Ostream& operator<<(Ostream& os, const chemistrySolution& l) {
    os << l.val;
    os << l.cellid;
    return os;
}
Istream& operator>>(Istream& is, chemistrySolution& l) {
    is >> l.val;
    is >> l.cellid;
    return is;
}

} // namespace Foam
