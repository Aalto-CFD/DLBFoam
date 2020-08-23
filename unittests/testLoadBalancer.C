#include "../third_party/catch.hpp"

#include "helpers.H"
#include "LoadBalancer.H"
#include "ChemistryProblem.H"


namespace Foam{




struct globalTest : public LoadBalancer {

public:
    using LoadBalancer::Operation;
    using LoadBalancer::getMin;
    using LoadBalancer::getMax;
    using LoadBalancer::getOperations;
    using LoadBalancer::timesToProblemCounts;
};


TEST_CASE("simpleBalancingMethod getMin()/getMax()"){

    DynamicList<ChemistryLoad> loads;
    loads.append(ChemistryLoad(1, 1.0));
    loads.append(ChemistryLoad(2, 2.0));
    loads.append(ChemistryLoad(3, 3.0));

    auto min = globalTest::getMin(loads);
    auto max = globalTest::getMax(loads);

    CHECK(min.rank == 1);
    CHECK(min.value == 1.0);
    CHECK(max.rank == 3);
    CHECK(max.value == 3);

}

TEST_CASE("LoadBalancer timesToProblemCounts"){


    size_t n_problems = 3;
    auto problems = create_random_problems(n_problems);

    set_cpu_times(problems, 1.0);
    std::vector<double> times = {1.1, 1.1};


    auto counts = globalTest::timesToProblemCounts(times, problems);

    CHECK(counts.size() == 3);
    CHECK(counts[0] == 1);
    CHECK(counts[1] == 1);
    CHECK(counts[2] == 1);


    set_cpu_times(problems, 1.0);
    times = {2.1, 0.9};

    counts = globalTest::timesToProblemCounts(times, problems);

    CHECK(counts.size() == 3);
    CHECK(counts[0] == 2);
    CHECK(counts[1] == 0);
    CHECK(counts[2] == 1);


}

TEST_CASE("LoadBalancer getOperations"){

    size_t n_tests = 50;

    for (size_t i = 0; i < n_tests; ++i){
        auto loads = create_random_load(i);

        for (size_t j = 0; j < i; ++j){
            REQUIRE_NOTHROW(globalTest::getOperations(loads, loads[i]));
        }
            
    }



}

/*


TEST_CASE("LoadBalancer updateState0()"){


    size_t n_problems = 100;
    size_t nprocs = Pstream::nProcs();
    
    
    double local_load_total;

    if (Pstream::master()){
        local_load_total = 1000.0;        
    }

    else {
        local_load_total = 10.0;
    }


    globalTest t;

    auto problems = get_problems_for_load(n_problems, local_load_total);

    auto myLoad   = t.computeLoad(problems);
    double global_mean = t.getMean(t.allGather(myLoad));


    


    t.updateState(problems);


    auto state = t.getState();

    if (Pstream::master()){
        CHECK(state.destinations.size() == nprocs);
        CHECK(state.sources.size() == 1);

        auto counts = state.nProblems;
        scalar remaining_load = counts[0] * local_load_total/n_problems;


        CHECK(std::abs(remaining_load - global_mean) < t.get_balancer_tolerance() * global_mean );

        size_t total = std::accumulate(counts.begin(), counts.end(), 0);
        CHECK(total == n_problems);

    }
    else {
        CHECK(state.destinations.size() == 1);
        CHECK(state.sources.size() == 2);
        CHECK(state.sources[1] == 0);
    }


}

TEST_CASE("LoadBalancer updateState1()"){


    size_t n_problems = 100;
    size_t nprocs = Pstream::nProcs();
    
    
    double local_load_total;

    if (Pstream::myProcNo() % 2){
        local_load_total = 1000.0;        
    }

    else {
        local_load_total = 10.0;
    }


    globalTest t;

    auto problems = get_problems_for_load(n_problems, local_load_total);

    auto myLoad   = t.computeLoad(problems);
    double global_mean = t.getMean(t.allGather(myLoad));


    


    t.updateState(problems);

    t.printState();

    auto state = t.getState();

    auto counts = state.nProblems;
    CHECK(std::accumulate(counts.begin(), counts.end(), 0) == n_problems);

    //senders
    if (Pstream::myProcNo() % 2){
        CHECK(state.destinations.size() > 1);
        CHECK(state.sources.size() == 1);


        auto counts = state.nProblems;
        scalar remaining_load = counts[0] * local_load_total/n_problems;


        CHECK(std::abs(remaining_load - global_mean) < t.get_balancer_tolerance() * global_mean );

    }
    //receivers
    else {
        CHECK(state.destinations.size() == 1);
        
        CHECK(state.sources.size() >= 1);
        //CHECK(state.sources[1] == 0);
    }

}


TEST_CASE("LoadBalancer updateState2()"){


    
    size_t n_problems = Pstream::myProcNo() * 41;
    double local_load_total = std::abs(sin(Pstream::myProcNo() * 0.2));

    auto problems = get_problems_for_load(n_problems, local_load_total);

    globalTest t;

    auto myLoad   = t.computeLoad(problems);
    auto all_loads = t.allGather(myLoad);

    double global_mean = t.getMean(all_loads);



    t.updateState(problems);


    auto state = t.getState();

    auto counts = state.nProblems;

    CHECK(std::accumulate(counts.begin(), counts.end(), 0) == n_problems);

    //sender
    if (counts.size() > 1) {
        scalar remaining_load = counts[0] * local_load_total/n_problems;
        CHECK(std::abs(remaining_load - global_mean) < t.get_balancer_tolerance() * global_mean );
    }
    

}

*/


}