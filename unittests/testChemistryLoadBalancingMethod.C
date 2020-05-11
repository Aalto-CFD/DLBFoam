#include "../third_party/catch.hpp"

#include "chemistryLoadBalancingMethod.H"
#include "chemistryProblem.H"


namespace Foam{

//create some arbitrary data
DynamicList<chemistryProblem> create_problems(int count){



    DynamicList<chemistryProblem> problems;
    
    for (int i = 0; i < count; ++i){
        //send data
        chemistryProblem p;
        
        scalarField s(10);
        s = 32.04; 
        p.c = s;
        p.Ti = 13.0;
        p.pi = 54.0;
        p.deltaTChem = 0.3;
        p.cellid = i;

        problems.append(p);

    }

    return problems;

}


class testableLoadBalancing : public chemistryLoadBalancingMethod{

public:


    //handy for testing
    void set_test_state(const std::vector<int>& sources, const std::vector<int>& destinations, const std::vector<int>& n_probs) {

        sendRecvInfo i;

        i.sources = sources;
        i.destinations = destinations;
        i.number_of_problems = n_probs;

        this->set_state(i);
  
    }


private:

    chemistryLoad get_my_load() const override{
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


    



    sendRecvInfo determine_state(const DynamicList<chemistryLoad>& loads) const override{
        
        
        sendRecvInfo i;
        /*
        i.sources{};
        i.destinations{1,2,3}
        i.number_of_problems{3,3,4};
        */
        return i;
        


    }
};

} //namespace Foam




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
    DynamicList<chemistryProblem> send_buffer = create_problems(10);
    DynamicList<chemistryProblem> recv_buffer;



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


}

TEST_CASE("chemistryLoadBalancingMethod get_send_buffer()"){

    

    using namespace Foam;

    testableLoadBalancing l;

    std::vector<int> sources;
    std::vector<int> destinations;
    std::vector<int> n_problems;

    

    int n_total_problems = 10;

    
    auto problems = create_problems(n_total_problems);
    
    sources = {};
    destinations = {1,2,3};
    n_problems = {3,3,2};

    l.set_test_state(sources, destinations, n_problems);
    
    auto buffer = l.get_send_buffer(problems);

    
    for (size_t i = 0; i < destinations.size(); ++i){
        CHECK(buffer[i].size() == n_problems[i]);
        
    }
    

    

}