#include "chemistryProblem.H"

namespace Foam{


Ostream& operator<<(Ostream& os, const chemistryProblem& p){

    /*
    scalarField c;
    scalar Ti;
    scalar pi;
    scalar deltaTChem;
    label cellid;
    */
    for (size_t i = 0; i < p.c.size(); ++i){
        os << p.c[i];
    }
    os << p.Ti;
    os << p.pi;
    os << p.deltaTChem;
    os << p.cellid;


    return os;
}


Istream& operator>>(Istream& is, chemistryProblem& p){

    for (size_t i = 0; i < p.c.size(); ++i){
        is >> p.c[i];
    }
    is >> p.Ti;
    is >> p.pi;
    is >> p.deltaTChem;
    is >> p.cellid;


    return is;



    return is;

}


}