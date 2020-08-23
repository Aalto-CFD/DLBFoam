#include "../third_party/catch.hpp"

#include <vector>
#include "ChemistryProblem.H"
#include "algorithms_aalto.H"
#include "helpers.H"



TEST_CASE("Algorithms count_while for vector"){

    using namespace Foam;

    std::vector<int> v = {1, 10, 3, 5};

    int sum;
    int max;
    auto sum_upto = [&](const int i) {
        sum += i;
        return sum <= max;
    };

    //forward iterators
    sum = 0;
    max = 4;
    CHECK(count_while(v.begin(), v.end(), sum_upto) == 1);
    
    sum = 0;
    max = 11;
    CHECK(count_while(v.begin(), v.end(), sum_upto) == 2);

    sum = 0;
    max = 19;
    CHECK(count_while(v.begin(), v.end(), sum_upto) == 4);

    sum = 0;
    max = 1000;
    CHECK(count_while(v.begin(), v.end(), sum_upto) == 4);

    //reverse iterators
    sum = 0;
    max = 5;
    CHECK(count_while(v.rbegin(), v.rend(), sum_upto) == 1);

    sum = 0;
    max = 7;
    CHECK(count_while(v.rbegin(), v.rend(), sum_upto) == 1);

    sum = 0;
    max = 18;
    CHECK(count_while(v.rbegin(), v.rend(), sum_upto) == 3);

    sum = 0;
    max = 1000;
    CHECK(count_while(v.rbegin(), v.rend(), sum_upto) == 4);

}


TEST_CASE("Algorithms count_while for list of problems (forward iterators)"){


    using namespace Foam;

    size_t n_problems = 4;
    auto v = create_random_problems(n_problems);

    v[0].cpuTime = 1.0;
    v[1].cpuTime = 10.0;
    v[2].cpuTime = 3.0;
    v[3].cpuTime = 5.0;
    


    double sum;
    double max;
    auto sum_upto = [&](const ChemistryProblem& p) {
        sum += p.cpuTime;
        return sum <= max;
    };

    //forward iterators
    sum = 0.0;
    max = 4.0;
    CHECK(count_while(v.begin(), v.end(), sum_upto) == 1);

    sum = 0.0;
    max = 11.0;
    CHECK(count_while(v.begin(), v.end(), sum_upto) == 2);

    sum = 0.0;
    max = 19.0;
    CHECK(count_while(v.begin(), v.end(), sum_upto) == 4);

    sum = 0.0;
    max = 1000.0;
    CHECK(count_while(v.begin(), v.end(), sum_upto) == 4);
}

//THIS TEST CASE CATCHES A BUG THAT PETTERI SPENT SOME TIME AT.
//THE BEHAVIOUR OF std::reverse_iterator and Foam::DynamicList reverse iterators appear to be different.
//IT IS REQUIRED TO CALL (and implement) make_reverse TO GET A REAL REVERSE ITERATOR
// THE BUG IS FIXED IN OF DEV AS OF 23.8.2020
TEST_CASE("Algorithms count_while for list of problems (reverse iterators)"){


    using namespace Foam;

    size_t n_problems = 4;
    auto v = create_random_problems(n_problems);

    v[0].cpuTime = 1.0;
    v[1].cpuTime = 10.0;
    v[2].cpuTime = 3.0;
    v[3].cpuTime = 5.0;
    


    double sum;
    double max;
    auto sum_upto = [&](const ChemistryProblem& p) {
        sum += p.cpuTime;
        return sum <= max;
    };
    
     
    //reverse iterators
    sum = 0;
    max = 5;
    CHECK(count_while(make_reverse(v.end()), make_reverse(v.begin()), sum_upto) == 1);

    sum = 0;
    max = 7;
    CHECK(count_while(make_reverse(v.end()), make_reverse(v.begin()), sum_upto) == 1);

    sum = 0;
    max = 18;
    CHECK(count_while(make_reverse(v.end()), make_reverse(v.begin()), sum_upto) == 3);

    sum = 0;
    max = 1000;
    CHECK(count_while(make_reverse(v.end()), make_reverse(v.begin()), sum_upto) == 4);

}




