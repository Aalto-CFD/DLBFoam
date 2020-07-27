#pragma once

#include <stdlib.h>     //srand, rand

#include "simpleBalancingMethod.H"
#include "chemistryProblem.H"


double random_double(double min, double max);

Foam::DynamicList<Foam::chemistryProblem> create_random_problems(size_t count = 10);


Foam::DynamicList<Foam::chemistryLoad> create_random_load(size_t n);

//create some arbitrary data
static Foam::DynamicList<Foam::chemistryProblem> create_problems2(int count){

    using namespace Foam;

    DynamicList<chemistryProblem> problems;
    
    for (int i = 0; i < count; ++i){
        chemistryProblem p;
        
        scalarField s(10);
        s = 32.04; 
        p.c = s;
        p.Ti = 13.0;
        p.pi = 54.0;
        p.deltaTChem = 0.3;
        p.cellid = i;
        p.cpuTime = 1.1;
        problems.append(p);

    }

    return problems;

}
