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
    using globalBalancingMethod::get_operation;
    using globalBalancingMethod::get_operations;
};

static void print_ops(const std::vector<globalTest::Operation>& ops) {

    for (const auto& op : ops) {
        Info << "From: " << op.from.rank << " To: " << op.to.rank << " Value: " << op.value << endl;
    }

}

TEST_CASE("globalBalancingMethod get_operation()"){

    chemistryLoad sender(1, 200.0);
    chemistryLoad receiver(2, 100.0);

    double mean = 150.0;



    auto op = globalTest::get_operation(sender, receiver, mean);

    CHECK(op.from.rank == 1);
    CHECK(op.to.rank == 2);
    CHECK(op.value == 50.0);

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



TEST_CASE("globalBalancingMethod get_operations()"){


    auto loads = create_random_load(50);
    auto my_load = loads[0];
    double mean = globalTest::get_mean(loads);

    Info << "mean " << mean << endl;

    //print_loads(loads);
    auto ops = globalTest::get_operations(loads, mean, my_load);
    print_loads(loads);
    print_ops(ops);

    /*
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
    */
}


}