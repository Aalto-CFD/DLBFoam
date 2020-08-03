#pragma once

#include <stdlib.h>     //srand, rand

#include "simpleBalancingMethod.H"
#include "chemistryProblem.H"


double random_double(double min, double max);

Foam::DynamicList<Foam::chemistryProblem> create_random_problems(size_t count = 10);

Foam::DynamicList<Foam::chemistryLoad> create_random_load(size_t n);


void set_cpu_times(Foam::DynamicList<Foam::chemistryProblem>& problems, double cpu_time);
void print_loads(const Foam::DynamicList<Foam::chemistryLoad>& loads);