#include "chemistryLoad.H"

namespace Foam{




Ostream& operator<<(Ostream& os, const chemistryLoad& l){

    os << l.rank;
    os << l.value;
    os << l.number_of_active_cells;

    return os;
    

}

Istream& operator>>(Istream& is, chemistryLoad& l){

    is >> l.rank;
    is >> l.value;
    is >> l.number_of_active_cells;

    return is;
}

}