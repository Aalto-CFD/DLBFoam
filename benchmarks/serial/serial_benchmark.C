#define PICOBENCH_IMPLEMENT
#include "../../third_party/picobench.hpp"





static void asdf(picobench::state& s)
{

    picobench::scope scope(s); // Constructor starts measurement. Destrucror stops it
    for (int i=0; i<s.iterations(); ++i){
        double asd = 43.433 - 4343434.44 * 4444;
    }

}

PICOBENCH(asdf);





int main(int argc, char* argv[]){


    //ode<aaltoChemistryModelBase<psiReactionThermo, gasHThermoPhysics>> asd(g_state.thermo);

    picobench::runner runner;
    // Disregard command-line for simplicity

    // Two sets of iterations
    runner.set_default_state_iterations({1});

    // One sample per benchmark because the huge numbers are expected to compensate
    // for external factors
    runner.set_default_samples(10);

    // Run the benchmarks with some seed which guarantees the same order every time
//    auto report = runner.run_benchmarks(123);

    return runner.run();

    //runner.run_benchmarks(123);

    // Output to some file
//    report.to_csv(std::ofstream("my.csv"));

    return 0;


//    reaction->correct();


    





}





