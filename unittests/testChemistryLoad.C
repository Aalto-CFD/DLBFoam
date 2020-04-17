#include "../third_party/catch.hpp"
#include "chemistryLoad.H"


TEST_CASE("chemistryLoad comparisons"){


    using namespace Foam;

    chemistryLoad l1(1, 1.0);
    chemistryLoad l2(1, 2.0);

    /*
    CHECK(!(l1 == l2));
    CHECK(l1 != l2);
    CHECK(l1 <= l2);
    CHECK(l2 >= l1);
    CHECK(l1 < l2);
    CHECK(l2 > l2);
    */    


}