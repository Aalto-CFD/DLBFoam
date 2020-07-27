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

        for (auto& problem : problems){
            problem.cpuTime = 3.0;
        }

        if (n_problems > 0){
            CHECK(problems[0].cpuTime == 3.0);
        }


        auto load = TESTABLE::compute_my_load(problems);

        CHECK(load.value == problems.size() * 3.0);
        CHECK(load.rank == Pstream::myProcNo());
    }


    

}

TEST_CASE("simpleBalancingMethod times_to_problem_counts"){

    auto problems = create_problems2(3);
    auto one = TESTABLE::times_to_problem_counts({1.1}, problems);
    CHECK(one.size() == 1);
    CHECK(one[0] == 1);



    auto two = TESTABLE::times_to_problem_counts({1.1, 1.1}, problems);
    CHECK(two.size() == 2);
    CHECK(two[0] == 1);
    CHECK(two[1] == 1);

    auto three = TESTABLE::times_to_problem_counts({1.1, 1.1, 1.1}, problems);
    CHECK(three.size() == 3);
    CHECK(three[0] == 1);
    CHECK(three[1] == 1);
    CHECK(three[2] == 1);

}


TEST_CASE("simpleBalancingMethod build_tree()"){
    

    for (size_t n_nodes = 0; n_nodes < 20; ++n_nodes){

        auto loads = create_random_load(n_nodes);

        auto root = TESTABLE::build_tree(loads);
        CHECK(loadTree::find(root, -1) != nullptr);

        
        for (size_t i = 0; i < n_nodes; ++i) {
            CHECK(loadTree::find(root, i) != nullptr);
        }

        CHECK(loadTree::find(root, n_nodes + 1) == nullptr);
    }


    




}

TEST_CASE("simpleBalancingMethod update_state()"){

    
    simpleBalancingMethod l;

    for (size_t n_problems = 0; n_problems < 30; ++n_problems){
        auto problems = create_problems2(n_problems);
        REQUIRE_NOTHROW(l.update_state(problems));

        //l.update_state(problems);
    }


}



}