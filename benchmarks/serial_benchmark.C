#define PICOBENCH_IMPLEMENT
#include "../third_party/picobench.hpp"


#include "create_state.H"

static ChemState<parallel::off> g_state;

struct InitialCondition1{

    static void set() {


        Random  randomizer1(541); 
        Random  randomizer2(44); 
        Random  randomizer3(10);

        forAll(g_state.rho, i) {
            g_state.rho[i] = 1.2 + randomizer1.sample01<scalar>() * 3;
            g_state.thermo.T()[i] = 500 + randomizer2.sample01<scalar>() * 100;
        }

        forAll(g_state.rho, celli) {

            for (label i = 0; i < g_state.Y.size(); i++) { 
                g_state.Y[i][celli] = randomizer3.sample01<scalar>();
                //c_[i] = Y_[i][celli]; 
            }

        }

    }

};



static void sanity_check(picobench::state& s)
{

    

    picobench::scope scope(s); // Constructor starts measurement. Destrucror stops it
    for (int i=0; i<s.iterations(); ++i){

        InitialCondition1::set();
        auto r1  = g_state.model1->solve(1E-3);
        InitialCondition1::set();
        auto r2  = g_state.model2->solve(1E-3);
        
        if (r1 != r2) {
            throw "solutions differ exiting";
        }



    }

}

PICOBENCH(sanity_check);

static void standard_chemistry_solve(picobench::state& s)
{

    picobench::scope scope(s); // Constructor starts measurement. Destrucror stops it
    for (int i=0; i<s.iterations(); ++i){

        InitialCondition1::set();
        g_state.model2->solve(1E-3);
    }

}

PICOBENCH(standard_chemistry_solve);

static void pyjac_chemistry_solve(picobench::state& s)
{


    picobench::scope scope(s); // Constructor starts measurement. Destrucror stops it
    for (int i=0; i<s.iterations(); ++i){
    
        InitialCondition1::set();
        g_state.model1->solve(1E-3);
    }
}

PICOBENCH(pyjac_chemistry_solve);





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





