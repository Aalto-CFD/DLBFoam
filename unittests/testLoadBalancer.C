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


TEST_CASE("LoadBalancer constructors"){

    REQUIRE_NOTHROW(LoadBalancer());


}


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

    CHECK(counts.size() == 2);
    CHECK(counts[0] == 1);
    CHECK(counts[1] == 1);


    set_cpu_times(problems, 1.0);
    times = {2.1, 0.9};

    counts = globalTest::timesToProblemCounts(times, problems);

    CHECK(counts.size() == 2);
    CHECK(counts[0] == 2);
    CHECK(counts[1] == 0);



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



}