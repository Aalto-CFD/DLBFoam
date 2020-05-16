#include "chemistryProblem.H"

namespace Foam {

Ostream& operator<<(Ostream& os, const chemistryProblem& p) {

    os << p.c;
    os << p.Ti;
    os << p.pi;
    os << p.rhoi;
    os << p.deltaTChem;
    os << p.deltaT;
    os << p.cellid;

    return os;
}

Istream& operator>>(Istream& is, chemistryProblem& p) {
    
    is >> p.c;
    is >> p.Ti;
    is >> p.pi;
    is >> p.rhoi;
    is >> p.deltaTChem;
    is >> p.deltaT;
    is >> p.cellid;

    return is;
}

} // namespace Foam