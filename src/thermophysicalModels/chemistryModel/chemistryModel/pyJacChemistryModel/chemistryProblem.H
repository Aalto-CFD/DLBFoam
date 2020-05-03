#pragma once

#include "volFields.H"

namespace Foam{



struct chemistryProblem{
  
    scalarField c;
    scalar Ti;
    scalar pi;
    scalar deltaTChem;
    label cellid;

    //TODO: implement!
    bool operator==(const chemistryProblem& rhs) const{
        return false;
    }

    bool operator !=(const chemistryProblem& rhs) const{
        return !(*this == rhs);
    }

};

//These are required for sending and receving using the OF MPI routines
Ostream& operator<<(Ostream& os, const chemistryProblem& l);
Istream& operator>>(Istream& is, chemistryProblem& l);




} //namespace Foam