#include <vector>

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"


///
#include "thermoPhysicsTypes.H"
#include "Random.H"
#include "loadBalancedChemistryModel.H"
#include "StandardChemistryModel.H"
#include "noChemistrySolver.H"
#include "ode.H"

#include "initial_conditions.H"
#include "result.H"
#include "runner.H"
#include "benchmarks.H"







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


    thermo.correct(); 

    std::vector<Result> results;

    /*
    results.push_back(
        Runner::run( BenchmarkSolveSingle( "LoadBalanced solve_single() light", thermo, false), 200 )
                    );


    results.push_back(
        Runner::run( BenchmarkSolveSingle( "LoadBalanced solve_single() heavy", thermo, true), 200 )
                    );

    */


    assign_heavy(p, rho, Y, thermo);


    results.push_back(
        Runner::run(
            BenchmarkSolve("Standard solve() heavy", ModelType::standard, thermo, HighMasterLoadIc()),
            10
        )
    );
    
    assign_heavy(p, rho, Y, thermo);

    results.push_back(
        Runner::run(
            BenchmarkSolve("Loadbalanced solve() heavy", ModelType::balanced, thermo, HighMasterLoadIc()),
            10
        )
    );

    assign_light(p, rho, Y, thermo);


    results.push_back(
        Runner::run(
            BenchmarkSolve("Standard solve() light", ModelType::standard, thermo, HighMasterLoadIc()),
            10
        )
    );
    
    assign_light(p, rho, Y, thermo);

    results.push_back(
        Runner::run(
            BenchmarkSolve("Loadbalanced solve() light", ModelType::balanced, thermo, HighMasterLoadIc()),
            10
        )
    );

    /*
    auto result1 = Runner::run(
        BenchmarkSolve( "Standard solve()",
                        ModelType::standard, 
                        thermo, 
                        HighMasterLoadIc()), 
                        10);
    Info << result1.to_string() << endl;



    auto result2 = Runner::run(
        BenchmarkSolve( "Loadbalanced solve()",
                        ModelType::balanced, 
                        thermo, 
                        HighMasterLoadIc()), 
                        10);

    Info << result2.to_string() << endl;
    */


    for (auto r : results) {
        Info << r.to_string() << endl;
    }

    //Runner r1(B1<ModelType::standard>(m1), 10);


    /*
    auto times = Benchmark(B1(m1), 10).times;

    for (auto t : times) {
        Info << t << endl;
    }
    */

    return 0;


}

