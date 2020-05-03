#include "../third_party/catch.hpp"

#include "chemistryLoadBalancingMethod.H"
#include "chemistryProblem.H"


namespace Foam{

class testableLoadBalancing : public chemistryLoadBalancingMethod{

private:

    chemistryLoad get_my_load() const {
        chemistryLoad load;    
        load.rank = Pstream::myProcNo();
        load.number_of_active_cells = load.rank * 2;
        load.value = double(3.0 * load.rank);
        
        return load;
    }

    WHATTODO determine_state(const DynamicList<chemistryLoad>& loads) const {
        return WHATTODO::e_SENDER;
    }
};

}




TEST_CASE("chemistryLoadBalancingMethod constructors"){


    using namespace Foam;

    REQUIRE_NOTHROW(testableLoadBalancing());

}


TEST_CASE("chemistryLoadBalancingMethod get_loads()"){


    using namespace Foam;

    testableLoadBalancing l;

    auto loads = l.get_loads();
    CHECK(Pstream::nProcs() == loads.size());

    CHECK(std::is_sorted(loads.begin(), loads.end()));

    for (const auto & load : loads){
        
        CHECK(load.value == double(3.0 * load.rank));
        CHECK(load.number_of_active_cells == load.rank * 2);
        
    }
    
    
}


TEST_CASE("chemistryLoadBalancingMethod send()"){


    using namespace Foam;

    testableLoadBalancing l;

    

    DynamicList<chemistryProblem> problems;

    for (size_t i = 0; i < 10; ++i){
        chemistryProblem p;
        
        problems.append(p);
    
    }

    REQUIRE_NOTHROW(l.send_recv(problems, 0, 1));

    //Problem p;



}