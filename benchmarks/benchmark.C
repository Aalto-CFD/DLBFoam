#include <vector>


#include "fvCFD.H"
#include "fluidThermoMomentumTransportModel.H"
#include "psiReactionThermophysicalTransportModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
///
//#include "thermoPhysicsTypes.H"
#include "Random.H"
#include "LoadBalancedChemistryModel.H"
#include "StandardChemistryModel.H"
#include "noChemistrySolver.H"
#include "ode.H"

#include "benchmark_info.H"
#include "initial_conditions.H"
#include "result.H"
#include "runner.H"
#include "benchmarks.H"



void sanity_check(volScalarField& p, volScalarField& rho, PtrList<volScalarField>& Y, psiReactionThermo& thermo) {

    BasicChemistryModel<psiReactionThermo>* model1 = new ode<StandardChemistryModel<psiReactionThermo, gasHThermoPhysics>>(thermo);
    BasicChemistryModel<psiReactionThermo>* model2 = new ode<LoadBalancedChemistryModel<psiReactionThermo, gasHThermoPhysics>>(thermo);

    set_master_heavy(rho, thermo);

    for (size_t i = 0; i < 2; ++i){
        model1->solve(detail::FLOWTIMESTEP);
        model2->solve(detail::FLOWTIMESTEP);
    }



    for (label i = 0; i < model1->nSpecie(); ++i){
        
        
        const scalarField& r1 = model1->RR(i);
        const scalarField& r2 = model2->RR(i);


        scalar sum = 0.0;
        forAll(r1, celli) {

            scalar diff = r1[celli] - r2[celli];
            sum += (diff * diff);
        }

        scalar norm =std::sqrt(sum);


        if (norm > 1E-10){
            Info << "Load balanced model and standard model give different answers! Be careful!" << endl;
        }

    }

    delete model1; delete model2;


}

void dump_results(const std::vector<Result>& results) {


    std::string fname = "results_" + std::to_string(Pstream::myProcNo()) + ".dat";

    ofstream outfile;
    outfile.open(fname);

    outfile << Result::get_header_csv() << std::endl;

    //Foam::fileName outputFile(fname);
    //Foam::OFstream os(outputFile);

    //os << Result::get_header() << endl;

    for (const auto& r : results) {
        outfile << r.to_csv() << std::endl;
    }

}



int main(int argc, char *argv[])
{

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"


    #include "compressibleCourantNo.H"
    #include "setDeltaT.H"


    //thermo.correct(); 

    std::vector<Result> results;

    ///
    ///@brief Sanity check that models give same answer
    ///
    ///
    sanity_check(p, rho, Y, thermo);
    
    ///
    ///@brief Benchmark the load balanced solver for light and heavy problems
    ///
    ///
    /* 
    results.push_back(
        Runner::run( BenchmarkSolveSingle( {"loadBalanced", "solveSingle()",  "simple", "all heavy"}, thermo, false), 50 )
                    );

    results.push_back(
        Runner::run( BenchmarkSolveSingle( {"loadBalanced", "solveSingle()",  "simple", "all light"}, thermo, true), 50 )
                    );

    */


    /////////////////////////
    ///@brief All heavy problems on all ranks 
    ///
    ///

    
    set_all_heavy(rho, thermo);
    
    results.push_back(
        Runner::run(
            BenchmarkSolve({"Standard", "solve()","none", "all heavy"}, ModelType::standard, thermo),
            10
        )
    );

    results.push_back(
        Runner::run(
            BenchmarkSolve({"loadBalanced", "solve()","simple", "all heavy"}, ModelType::balanced, thermo),
            10
        )
    );

    ///
    ///@brief All light problems on all ranks
    ///
    ///

    set_all_light(rho, thermo);

    results.push_back(
        Runner::run(
            BenchmarkSolve({"Standard", "solve()","none", "all light"}, ModelType::standard, thermo),
            10
        )
    );
    

    results.push_back(
        Runner::run(
            BenchmarkSolve({"loadBalanced", "solve()","simple", "all light"}, ModelType::balanced, thermo),
            10
        )
    );


    ///
    ///@brief Set master to high value, others to light value
    ///
    ///
    set_master_heavy(rho, thermo);
    

    results.push_back(
        Runner::run(
            BenchmarkSolve({"Standard", "solve()","none", "master heavy"}, ModelType::standard, thermo),
            10
        )
    );
    

    results.push_back(
        Runner::run(
            BenchmarkSolve({"loadBalanced", "solve()", "simple", "master heavy"}, ModelType::balanced, thermo),
            10
        )
    );


    ///
    ///@brief Every second rank is heavy
    ///
    ///
    set_every_n_heavy(rho, thermo, 2);
    
    results.push_back(
        Runner::run(
            BenchmarkSolve({"Standard", "solve()","none", "every 2nd heavy"}, ModelType::standard, thermo),
            10
        )
    );
    

    results.push_back(
        Runner::run(
            BenchmarkSolve({"loadBalanced", "solve()","simple", "every 2nd heavy"}, ModelType::balanced, thermo),
            10
        )
    );


    ///
    ///@brief Every fourth rank is heavy
    ///
    ///
    set_every_n_heavy(rho, thermo, 4);
    
    results.push_back(
        Runner::run(
            BenchmarkSolve({"Standard", "solve()","none", "every 4th heavy"}, ModelType::standard, thermo),
            10
        )
    );
    

    results.push_back(
        Runner::run(
            BenchmarkSolve({"loadBalanced", "solve()","simple", "every 4th heavy"}, ModelType::balanced, thermo),
            10
        )
    );


    ///
    ///@brief Every tenth rank is heavy
    ///
    ///
    set_every_n_heavy(rho, thermo, 10);
    
    results.push_back(
        Runner::run(
            BenchmarkSolve({"Standard", "solve()","none", "every 10th heavy"}, ModelType::standard, thermo),
            10
        )
    );
    

    results.push_back(
        Runner::run(
            BenchmarkSolve({"loadBalanced", "solve()","simple", "every 10th heavy"}, ModelType::balanced, thermo),
            10
        )
    );

    ///
    ///@brief Every third rank is randomly heavy
    ///
    ///
    set_every_n_random_heavy(rho, thermo, 4);
    
    results.push_back(
        Runner::run(
            BenchmarkSolve({"Standard", "solve()","none", "every 3rd rnd heavy"}, ModelType::standard, thermo),
            10
        )
    );
    

    results.push_back(
        Runner::run(
            BenchmarkSolve({"loadBalanced", "solve()","simple", "every 3rd rnd heavy"}, ModelType::balanced, thermo),
            10
        )
    );


    ///
    ///@brief All random
    ///
    ///
    set_all_random(rho, thermo);
    
    results.push_back(
        Runner::run(
            BenchmarkSolve({"Standard", "solve()","none", "all random"}, ModelType::standard, thermo),
            10
        )
    );
    

    results.push_back(
        Runner::run(
            BenchmarkSolve({"loadBalanced", "solve()","simple", "all random"}, ModelType::balanced, thermo),
            10
        )
    );

    dump_results(results);

    

    Info << Result::get_header() << endl;

    for (auto r : results) {
        
        Info << r << endl;
    }

    return 0;


}

