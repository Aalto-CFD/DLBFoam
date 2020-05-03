#pragma once

#include "volFields.H"

namespace Foam{

struct chemistryProblem{
  
    scalarField c;
    scalar Ti;
    scalar pi;
    scalar deltaTChem;
    label cellid;
};


} //namespace Foam