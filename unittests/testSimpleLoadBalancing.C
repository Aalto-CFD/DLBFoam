#include "../third_party/catch.hpp"

#include "helpers.H"
#include "simpleBalancingMethod.H"
#include "chemistryProblem.H"


namespace Foam{




struct TESTABLE : public simpleBalancingMethod {

public:
    using simpleBalancingMethod::compute_my_load;
    using simpleBalancingMethod::build_tree;
    using simpleBalancingMethod::times_to_problem_counts;
};




TEST_CASE("simpleBalancingMethod compute_my_load()"){

    TESTABLE l;
    
    for (size_t n_problems = 0; n_problems < 30; ++n_problems){

        auto problems = create_random_problems(n_problems);

        set_cpu_times(problems, 3.0);

        if (n_problems > 0){
            CHECK(problems[0].cpuTime == 3.0);
        }


        auto load = TESTABLE::compute_my_load(problems);

        CHECK(load.value == problems.size() * 3.0);
        CHECK(load.rank == Pstream::myProcNo());
    }

}


TEST_CASE("simpleBalancingMethod times_to_problem_counts1"){


    size_t n_problems = 10;
    auto problems = create_random_problems(n_problems);

    set_cpu_times(problems, 1.0);

    auto counts1 = TESTABLE::times_to_problem_counts({5.1, 1.1, 1.1, 1.1, 1.1, 1.1}, problems);
    CHECK(counts1.size() == 6);
    CHECK(std::accumulate(counts1.begin(), counts1.end(), 0) == n_problems);

    auto counts2 = TESTABLE::times_to_problem_counts({2.0, 3.0, 5.0}, problems);
    CHECK(counts2.size() == 3);
    CHECK(std::accumulate(counts2.begin(), counts2.end(), 0) == n_problems);


    auto counts3 = TESTABLE::times_to_problem_counts({10.0}, problems);
    CHECK(counts3.size() == 1);
    CHECK(std::accumulate(counts3.begin(), counts3.end(), 0) == n_problems);


    //auto counts4 = TESTABLE::times_to_problem_counts({1.1}, problems);
    //CHECK(counts4.size() == 1);
    //CHECK(std::accumulate(counts4.begin(), counts4.end(), 0) == 1);

}


TEST_CASE("simpleBalancingMethod times_to_problem_counts2"){


    size_t n_problems = 3;
    auto problems = create_random_problems(n_problems);

    std::vector<double> times;
    times.push_back(problems[0].cpuTime);
    times.push_back(problems[1].cpuTime + problems[2].cpuTime);

    auto counts = TESTABLE::times_to_problem_counts(times, problems);

    CHECK(counts.size() == 2);
    CHECK(std::accumulate(counts.begin(), counts.end(), 0) == n_problems);

}

TEST_CASE("simpleBalancingMethod build_tree()"){
   /* 

    for (size_t n_nodes = 0; n_nodes < 20; ++n_nodes){

        auto loads = create_random_load(n_nodes);

        auto root = TESTABLE::build_tree(loads);
        CHECK(loadTree::find(root, -1) != nullptr);

        
        for (size_t i = 0; i < n_nodes; ++i) {
            CHECK(loadTree::find(root, i) != nullptr);
        }

        CHECK(loadTree::find(root, n_nodes + 1) == nullptr);

        //loadTree::print(root);
    }
    */

    




}
TEST_CASE("simpleBalancingMethod update_state()"){

    TESTABLE t;

    size_t n_problems = 500;

    auto problems = create_random_problems(n_problems);

    for (size_t i = 0; i < problems.size(); ++i) {
        problems[i].cpuTime = (Pstream::myProcNo() + 1) * 5 * (i+1); // * random_double(1, 100);
    }

    auto load = TESTABLE::compute_my_load(problems);
    CHECK(load.value > 0);

    REQUIRE_NOTHROW(t.update_state(problems));


}


} //namespace Foam