#include "../third_party/catch.hpp"

#include "simpleBalancingMethod.H"
#include "chemistryProblem.H"


namespace Foam{




struct TESTABLE : public simpleBalancingMethod {

public:
    using simpleBalancingMethod::compute_my_load;
    using simpleBalancingMethod::build_tree;
    using simpleBalancingMethod::times_to_problem_counts;
};


static Foam::DynamicList<Foam::chemistryLoad> create_random_load(size_t n){
    
    using namespace Foam;

    DynamicList<chemistryLoad> ret(n, chemistryLoad());

    for( size_t i = 0; i < n; ++i ) {
        ret[i].value = double(rand() % 100) + 1;
        ret[i].rank = i;
    }
    return ret;
}


//create some arbitrary data
DynamicList<chemistryProblem> create_problems2(int count){



    DynamicList<chemistryProblem> problems;
    
    for (int i = 0; i < count; ++i){
        chemistryProblem p;
        
        scalarField s(10);
        s = 32.04; 
        p.c = s;
        p.Ti = 13.0;
        p.pi = 54.0;
        p.deltaTChem = 0.3;
        p.cellid = i;
        p.cpuTime = 1.1;
        problems.append(p);

    }

    return problems;

}




TEST_CASE("simpleBalancingMethod compute_my_load()"){


    
    TESTABLE l;

    auto problems = create_problems2(3 + 1*Pstream::myProcNo());
    CHECK(problems[0].cpuTime == 1.1);

    auto load = TESTABLE::compute_my_load(problems);

    CHECK(load.value == problems.size() * 1.1);
    CHECK(load.rank == Pstream::myProcNo());

}

TEST_CASE("simpleBalancingMethod build_tree()"){
    
    size_t n_nodes = 20;

    auto loads = create_random_load(n_nodes);

    auto root = TESTABLE::build_tree(loads);
    
    CHECK(loadTree::find(root, -1) != nullptr);

    for (size_t i = 0; i < n_nodes; ++i) {
        CHECK(loadTree::find(root, i) != nullptr);
    }

    CHECK(loadTree::find(root, 21) == nullptr);

    //loadTree::print(root, n_nodes);

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


TEST_CASE("simpleBalancingMethod update_state()"){

    
    simpleBalancingMethod l;

    auto problems = create_problems2(10);
    l.update_state(problems);

    /*
    auto buffer_s = l.balance(problems);

    Pout << buffer_s.size() << endl;

    //CHECK(buffer_s.size() != 0);

    */
}



}