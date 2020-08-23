#include "helpers.H"

double random_double(double min, double max) {
    double f = double(rand()) / RAND_MAX;
    return min + f * (max - min);
}

Foam::DynamicList<Foam::ChemistryProblem> create_random_problems(size_t count) {

    using namespace Foam;
    DynamicList<ChemistryProblem> problems;

    for (size_t i = 0; i < count; ++i) {
        ChemistryProblem p;

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

Foam::DynamicList<Foam::ChemistryLoad> create_random_load(size_t count)
{
    using namespace Foam;

    DynamicList<ChemistryLoad> ret(count, ChemistryLoad());

    for (size_t i = 0; i < count; ++i) {
        ret[i].value = random_double(0., 300.0);
        ret[i].rank  = i;
    }
    return ret;
}

Foam::DynamicList<Foam::ChemistryProblem> getProblems_for_load(size_t n_problems, double total_load){

    double problem_value = total_load / n_problems;

    auto problems = create_random_problems(n_problems);


    for (auto& problem : problems){
        problem.cpuTime = problem_value;
    }
    return problems;

}

void set_cpu_times(Foam::DynamicList<Foam::ChemistryProblem>& problems, double cpu_time){

    for (auto& problem : problems){
        problem.cpuTime = cpu_time;
    }

}

void print_loads(const Foam::DynamicList<Foam::ChemistryLoad>& loads) {

    for (const auto& load : loads) {

        Foam::Info << "Rank: " << load.rank << " Value: " << load.value << Foam::endl;

    }

}
