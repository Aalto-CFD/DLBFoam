#include "../third_party/catch.hpp"

#include "chemistryLoadBalancingMethod.H"
#include "chemistryProblem.H"


namespace Foam{

class testableLoadBalancing : public chemistryLoadBalancingMethod{

private:

    chemistryLoad get_my_load() const {
        chemistryLoad load;    
        load.rank = Pstream::myProcNo();

        if (load.rank % 2 == 0){
            load.number_of_active_cells = load.rank * 2 + 1;
            load.value = double(3.0 * load.rank + 1);
        }

        else{
            load.number_of_active_cells = load.rank * 4 + 1;
            load.value = double(3.65 * load.rank + 1);
        }
        
        
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
        

        if (load.rank % 2 == 0){

            CHECK(load.value == double(3.0 * load.rank + 1));
            CHECK(load.number_of_active_cells == load.rank * 2 + 1);
        
        }

        else{
            CHECK(load.value == double(3.65 * load.rank + 1));
            CHECK(load.number_of_active_cells == load.rank * 4 + 1);
        
        }

    }

    
}



TEST_CASE("chemistryLoadBalancingMethod Isend_recv()"){


    using namespace Foam;

    testableLoadBalancing l;
    
    int source = 0;
    int destination = 1;
    DynamicList<chemistryProblem> send_buffer;
    DynamicList<chemistryProblem> recv_buffer;

    //create some arbitrary data
    for (size_t i = 0; i < 10; ++i){
        //send data
        chemistryProblem p;
        
        scalarField s(10);
        s = 32.04; 
        p.c = s;
        p.Ti = 13.0;
        p.pi = 54.0;
        p.deltaTChem = 0.3;
        p.cellid = i;

        send_buffer.append(p);

        //receive data, set to zero here so that checking can be done
        chemistryProblem p2;
        s = 0.0; 
        p2.c = s;
        p2.Ti = 0.0;
        p2.pi = 0.0;
        p2.deltaTChem = 0.0;
        p2.cellid = 0;

        recv_buffer.append(p2);
    }


    l.Isend_recv(send_buffer, recv_buffer, source, destination);

    if (Pstream::myProcNo() == destination){
        for (size_t i = 0; i < 10; ++i){
            auto p = recv_buffer[i];
            CHECK(p.c[0] == 32.04); //note only the first one checked
            CHECK(p.Ti == 13.0);
            CHECK(p.pi == 54.0);
            CHECK(p.deltaTChem == 0.3);
            CHECK(p.cellid == i);
        }
    }


    //REQUIRE_NOTHROW(l.send_recv(problems, 0, 1));

    //Problem p;



}