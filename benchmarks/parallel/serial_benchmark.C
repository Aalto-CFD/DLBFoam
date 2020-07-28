#define PICOBENCH_IMPLEMENT
#include "../../third_party/picobench.hpp"

#include <iostream>
#include <fstream> 

#include "create_state.H"
#include "initial_conditions.H"



static void sanity_check(picobench::state& s) {


    static ChemState<parallel::off, chemistryModelName::standard> g_state1;
    static ChemState<parallel::off, chemistryModelName::loadBalanced> g_state2;


    picobench::scope scope(s);

    auto r1 = g_state1.model->solve(1E-3);
    auto r2 = g_state2.model->solve(1E-3);

    Info << "Standard model solution: " << r1 << endl;
    Info << "Loadbalanced model solution: " << r2 << endl;


}

PICOBENCH(sanity_check);



static void standard_chemistry_solve(picobench::state& s)
{


    static ChemState<parallel::off, chemistryModelName::standard> g_state;
    random_rho_T_Y_initial_condition::set(g_state);

    picobench::scope scope(s); // Constructor starts measurement. Destrucror stops it
    for (int i=0; i<s.iterations(); ++i){

        auto r = g_state.model->solve(1E-3);
    }

}

PICOBENCH(standard_chemistry_solve);

static void loadbalanced_chemistry_solve(picobench::state& s)
{


    static ChemState<parallel::off, chemistryModelName::loadBalanced> g_state;
    random_rho_T_Y_initial_condition::set(g_state);

    picobench::scope scope(s); // Constructor starts measurement. Destrucror stops it
    for (int i=0; i<s.iterations(); ++i){

        auto r = g_state.model->solve(1E-3);
    }

}

PICOBENCH(loadbalanced_chemistry_solve);




int main(int argc, char* argv[]){



    picobench::runner runner;
    runner.parse_cmd_line(argc, argv);
    // Disregard command-line for simplicity

    // Two sets of iterations
    runner.set_default_state_iterations({1});

    // One sample per benchmark because the huge numbers are expected to compensate
    // for external factors
    runner.set_default_samples(10);

    runner.run();

    auto report = runner.generate_report();


    ofstream file;
    file.open ("benchmark_results.csv");
    report.to_csv(file); // Otputs in csv format. Most detailed
    //runner.run_benchmarks(123);

    // Output to some file
//    report.to_csv(std::ofstream("my.csv"));

    return 0;


//    reaction->correct();


    





}





