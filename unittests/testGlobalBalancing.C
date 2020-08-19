#include "../third_party/catch.hpp"

#include "helpers.H"
#include "globalBalancingMethod.H"
#include "chemistryProblem.H"


namespace Foam{




struct globalTest : public globalBalancingMethod {

public:
    using globalBalancingMethod::Operation;
    using globalBalancingMethod::get_min;
    using globalBalancingMethod::get_max;
    using globalBalancingMethod::get_operations;
    using globalBalancingMethod::times_to_problem_counts;
};

static void print_ops(const std::vector<globalTest::Operation>& ops) {

    for (const auto& op : ops) {
        Info << "From: " << op.from.rank << " To: " << op.to.rank << " Value: " << op.value << endl;
    }

}

TEST_CASE("simpleBalancingMethod get_min()/get_max()"){

    DynamicList<chemistryLoad> loads;
    loads.append(chemistryLoad(1, 1.0));
    loads.append(chemistryLoad(2, 2.0));
    loads.append(chemistryLoad(3, 3.0));

    auto min = globalTest::get_min(loads);
    auto max = globalTest::get_max(loads);

    CHECK(min.rank == 1);
    CHECK(min.value == 1.0);
    CHECK(max.rank == 3);
    CHECK(max.value == 3);

}
/*
TEST_CASE("globalBalancingMethod times_to_problem_counts1"){


    size_t n_problems = 10;
    auto problems = create_random_problems(n_problems);

    set_cpu_times(problems, 1.0);

    auto counts1 = globalTest::times_to_problem_counts({5.1, 1.1, 1.1, 1.1, 1.1, 1.1}, problems);
    CHECK(counts1.size() == 6);
    CHECK(std::accumulate(counts1.begin(), counts1.end(), 0) == n_problems);

    auto counts2 = globalTest::times_to_problem_counts({2.0, 3.0, 5.0}, problems);
    CHECK(counts2.size() == 3);
    CHECK(std::accumulate(counts2.begin(), counts2.end(), 0) == n_problems);


    auto counts3 = globalTest::times_to_problem_counts({10.0}, problems);
    CHECK(counts3.size() == 1);
    CHECK(std::accumulate(counts3.begin(), counts3.end(), 0) == n_problems);


    //auto counts4 = TESTABLE::times_to_problem_counts({1.1}, problems);
    //CHECK(counts4.size() == 1);
    //CHECK(std::accumulate(counts4.begin(), counts4.end(), 0) == 1);

}
*/

TEST_CASE("globalBalancingMethod times_to_problem_counts"){


    size_t n_problems = 3;
    auto problems = create_random_problems(n_problems);

    set_cpu_times(problems, 1.0);
    std::vector<double> times = {1.1, 1.1};


    auto counts = globalTest::times_to_problem_counts(times, problems);

    CHECK(counts.size() == 3);
    CHECK(counts[0] == 1);
    CHECK(counts[1] == 1);
    CHECK(counts[2] == 1);


    set_cpu_times(problems, 1.0);
    times = {2.1, 0.9};

    counts = globalTest::times_to_problem_counts(times, problems);

    CHECK(counts.size() == 3);
    CHECK(counts[0] == 2);
    CHECK(counts[1] == 0);
    CHECK(counts[2] == 1);


}

/*


TEST_CASE("globalBalancingMethod update_state0()"){


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

    auto my_load   = t.compute_my_load(problems);
    double global_mean = t.get_mean(t.all_gather(my_load));


    


    t.update_state(problems);


    auto state = t.get_state();

    if (Pstream::master()){
        CHECK(state.destinations.size() == nprocs);
        CHECK(state.sources.size() == 1);

        auto counts = state.number_of_problems;
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

TEST_CASE("globalBalancingMethod update_state1()"){


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

    auto my_load   = t.compute_my_load(problems);
    double global_mean = t.get_mean(t.all_gather(my_load));


    


    t.update_state(problems);

    t.print_state();

    auto state = t.get_state();

    auto counts = state.number_of_problems;
    CHECK(std::accumulate(counts.begin(), counts.end(), 0) == n_problems);

    //senders
    if (Pstream::myProcNo() % 2){
        CHECK(state.destinations.size() > 1);
        CHECK(state.sources.size() == 1);


        auto counts = state.number_of_problems;
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


TEST_CASE("globalBalancingMethod update_state2()"){


    
    size_t n_problems = Pstream::myProcNo() * 41;
    double local_load_total = std::abs(sin(Pstream::myProcNo() * 0.2));

    auto problems = get_problems_for_load(n_problems, local_load_total);

    globalTest t;

    auto my_load   = t.compute_my_load(problems);
    auto all_loads = t.all_gather(my_load);

    double global_mean = t.get_mean(all_loads);



    t.update_state(problems);


    auto state = t.get_state();

    auto counts = state.number_of_problems;

    CHECK(std::accumulate(counts.begin(), counts.end(), 0) == n_problems);

    //sender
    if (counts.size() > 1) {
        scalar remaining_load = counts[0] * local_load_total/n_problems;
        CHECK(std::abs(remaining_load - global_mean) < t.get_balancer_tolerance() * global_mean );
    }
    

}

*/


}