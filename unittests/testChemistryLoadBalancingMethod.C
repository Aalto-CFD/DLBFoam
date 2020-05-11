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


    void update() {

        auto loads = get_loads();
        auto state = determine_state(loads);
        this->set_state(state);
    }

private:

    chemistryLoad get_my_load() const override{
        chemistryLoad load;    
        load.rank = Pstream::myProcNo();

        if (load.rank % 2 == 0){
            load.number_of_active_cells = load.rank * 2 + 50;
            load.value = double(3.0 * load.rank + 50);
        }

        else{
            load.number_of_active_cells = load.rank * 4 + 50;
            load.value = double(3.65 * load.rank + 50);
        }
        
        
        return load;
    }



    sendRecvInfo determine_state(const DynamicList<chemistryLoad>& loads) const override{
        
        
        sendRecvInfo i;


        if (Pstream::myProcNo() == Pstream::master()){

            for (const auto& load : loads){
                if (load.rank != Pstream::myProcNo()){

                    i.sources.emplace_back(load.rank);
                    i.number_of_problems.emplace_back(10);
                }
            }
            i.destinations = {};
        }

        else {

            /*for (const auto& load : loads){
                if (load.rank != Pstream::myProcNo()){

                    i.number_of_problems.emplace_back(10);
                }
            }*/
            i.number_of_problems = {10};
            i.destinations = {Pstream::master()};
            i.sources = {};
        }



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

            CHECK(load.value == double(3.0 * load.rank + 50));
            CHECK(load.number_of_active_cells == load.rank * 2 + 50);
        
        }

        else{
            CHECK(load.value == double(3.65 * load.rank + 50));
            CHECK(load.number_of_active_cells == load.rank * 4 + 50);
        
        }

    }

    
}



TEST_CASE("chemistryLoadBalancingMethod send_recv()"){

    
    using namespace Foam;

    testableLoadBalancing l;
    std::vector<int> sources = {};
    std::vector<int> destinations = {};
    std::vector<int> counts = {};

    if (Pstream::myProcNo() == 0){
        destinations = {1};
        counts = {10};
    }
    
    else if (Pstream::myProcNo() == 1){
        sources = {0};
        counts = {10};
    }

  

    DynamicList<DynamicList<chemistryProblem>> send_buffer;
    send_buffer.setSize(1);
    send_buffer[0] = create_problems(10);

    auto recv_buffer = l.send_recv<chemistryProblem, Pstream::commsTypes::nonBlocking>(
        send_buffer, sources, destinations);

    if (Pstream::myProcNo() == 1) {
        for (size_t i = 0; i < 10; ++i) {
            auto p = recv_buffer[0][i];
            CHECK(p.c[0] == 32.04); // note only the first one checked
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
    
    //cellids 0,1,2 "remain" in problems the first problem to send contains
    //cell id 2 
    CHECK(buffer[0][0].cellid == 2);
    CHECK(buffer[1][0].cellid == 2 + n_problems[0]);
    CHECK(buffer[2][0].cellid == 2 + n_problems[0] + n_problems[1]);


    n_problems = {4,4,4};
    l.set_test_state(sources, destinations, n_problems);
    
    REQUIRE_THROWS(l.get_send_buffer(problems));
    

}


/*
TEST_CASE("chemistryLoadBalancingMethod get_problems()"){

    

    using namespace Foam;

    

    testableLoadBalancing l;

    l.update();

    auto problems = create_problems(50);


    //l.set_state(l.determine_state(loads));

    auto buffer = l.get_problems(problems);


}
*/