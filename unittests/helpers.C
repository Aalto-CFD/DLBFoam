#include "helpers.H"

double random_double(double min, double max) {
    double f = double(rand()) / RAND_MAX;
    return min + f * (max - min);
}

Foam::DynamicList<Foam::chemistryProblem> create_random_problems(size_t count) {

    using namespace Foam;
    DynamicList<chemistryProblem> problems;

    for (size_t i = 0; i < count; ++i) {
        chemistryProblem p;

        scalarField s(10);
        p.c          = random_double(0.0, 1.0);
        p.Ti         = random_double(0.0, 1.0);
        p.pi         = random_double(0.0, 1.0);
        p.deltaTChem = random_double(0.0, 1.0);
        p.cellid     = i;
        p.cpuTime    = random_double(0.0, 1.0);
        problems.append(p);
    }

    return problems;
}

Foam::DynamicList<Foam::chemistryLoad> create_random_load(size_t count)
{
    using namespace Foam;

    DynamicList<chemistryLoad> ret(count, chemistryLoad());

    for (size_t i = 0; i < count; ++i) {
        ret[i].value = random_double(0., 300.0);
        ret[i].rank  = i;
    }
    return ret;
}


void set_cpu_times(Foam::DynamicList<Foam::chemistryProblem>& problems, double cpu_time){

    for (auto& problem : problems){
        problem.cpuTime = cpu_time;
    }

}

void print_loads(const Foam::DynamicList<Foam::chemistryLoad>& loads) {

    for (const auto& load : loads) {

        Foam::Info << "Rank: " << load.rank << " Value: " << load.value << Foam::endl;

    }

}