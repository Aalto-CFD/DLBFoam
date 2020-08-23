#pragma once

#include <stdlib.h>     //srand, rand

#include "ChemistryLoad.H"
#include "ChemistryProblem.H"


double random_double(double min, double max);

Foam::DynamicList<Foam::ChemistryProblem> create_random_problems(size_t count = 10);

Foam::DynamicList<Foam::ChemistryLoad> create_random_load(size_t n);
Foam::DynamicList<Foam::ChemistryProblem> getProblems_for_load(size_t n_problems, double total_load);

void set_cpu_times(Foam::DynamicList<Foam::ChemistryProblem>& problems, double cpu_time);
void print_loads(const Foam::DynamicList<Foam::ChemistryLoad>& loads);
