#include "chemistryLoad.H"

namespace Foam{




Ostream& operator<<(Ostream& os, const chemistryLoad& l){

    os << l.rank;
    os << l.value;

    return os;
    

}

Istream& operator>>(Istream& is, chemistryLoad& l){

    is >> l.rank;
    is >> l.value;

    return is;
}

}