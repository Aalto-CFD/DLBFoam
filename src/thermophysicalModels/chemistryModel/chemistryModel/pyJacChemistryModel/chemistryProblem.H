#pragma once

#include "volFields.H"

namespace Foam{



struct chemistryProblem{
  
    scalarField c;
    scalar Ti;
    scalar pi;
    scalar deltaTChem;
    label cellid;

    bool operator==(const chemistryProblem& rhs) const{
        return false;
    }

    bool operator !=(const chemistryProblem& rhs) const{
        return !(*this == rhs);
    }

};

Ostream& operator<<(Ostream& os, const chemistryProblem& p);




} //namespace Foam