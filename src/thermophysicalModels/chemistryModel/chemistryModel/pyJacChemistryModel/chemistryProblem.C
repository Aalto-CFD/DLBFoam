#include "chemistryProblem.H"

namespace Foam {

Ostream& operator<<(Ostream& os, const chemistryProblem& p) {

    os << p.c;
    os << p.Ti;
    os << p.pi;
    os << p.deltaTChem;
    os << p.cellid;

    return os;
}

Istream& operator>>(Istream& is, chemistryProblem& p) {
    
    is >> p.c;
    is >> p.Ti;
    is >> p.pi;
    is >> p.deltaTChem;
    is >> p.cellid;

    return is;
}

} // namespace Foam