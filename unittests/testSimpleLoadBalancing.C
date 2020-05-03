#include "../third_party/catch.hpp"

#include "simpleLoadBalancing.H"
#include "chemistryProblem.H"

TEST_CASE("simpleLoadBalancing constructors"){


    using namespace Foam;

    simpleLoadBalancing l;
    

}


TEST_CASE("simpleLoadBalancing get_loads()"){


    using namespace Foam;

    simpleLoadBalancing l;

    auto loads = l.get_loads();
    CHECK(Pstream::nProcs() == loads.size());

    CHECK(std::is_sorted(loads.begin(), loads.end()));

    /*
    for (auto load : loads){
        Info << load.rank << " " << load.value << endl;
    }
    */
}


TEST_CASE("simpleLoadBalancing send()"){


    using namespace Foam;

    simpleLoadBalancing l;

    //Problem p;



}