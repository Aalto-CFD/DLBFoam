#include "../third_party/catch.hpp"

#include "LoadBalancerBase.H"
#include "ChemistryProblem.H"


namespace Foam{

//create some arbitrary data
DynamicList<ChemistryProblem> create_problems(int count){

    DynamicList<ChemistryProblem> problems;
    
    for (int i = 0; i < count; ++i){
        //send data
        ChemistryProblem p;
        
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


class testableLoadBalancing : public LoadBalancerBase{

public:


    //handy for testing
    void set_test_state(const std::vector<int>& sources, const std::vector<int>& destinations, const std::vector<int>& n_probs) {

        BalancerState i;

        i.sources = sources;
        i.destinations = destinations;
        i.nProblems = n_probs;

        this->setState(i);
  
    }
    
    void updateState(const DynamicList<ChemistryProblem>& problems) {
        /*empty*/
    }

    template<class ET, Pstream::commsTypes CT>
    static buffer_t<ET> test_sendRecv(const buffer_t<ET>& buffer, std::vector<int> sources, std::vector<int> dests){
        return LoadBalancerBase::sendRecv<ET, CT>(buffer, sources, dests);
    }
   

   /* void update() {

        auto loads = get_loads();
        auto state = determine_state(loads);
        this->setState(state);
    }

    */
};

} //namespace Foam




TEST_CASE("LoadBalancerBase constructors"){


    using namespace Foam;

    REQUIRE_NOTHROW(testableLoadBalancing());

}

TEST_CASE("LoadBalancerBase allGather()"){



    using namespace Foam;

    double myData = 2.5 * Pstream::myProcNo();

    auto data = LoadBalancerBase::allGather(myData);

    CHECK(data[Pstream::myProcNo()] == 2.5 * Pstream::myProcNo());


    for (int i = 0; i < data.size(); ++i){

        CHECK(data[i] == double(i * 2.5));

    }

    CHECK(data.size() == Pstream::nProcs());
    


}



TEST_CASE("LoadBalancerBase partition()"){


    using namespace Foam;

    DynamicList<int> arr(10, int());

    for (int i = 0; i < 10; ++i){
        arr[i] = i;
    }
    std::vector<int> subsizes{2, 2, 1, 5};
 
    auto r = LoadBalancerBase::partition(arr, subsizes);

    CHECK(r[0][0] == 0);
    CHECK(r[0][1] == 1);
    CHECK(r[1][0] == 2);
    CHECK(r[3][4] == 9);


    auto problems = create_problems(10);

    auto split = LoadBalancerBase::partition(problems, subsizes);

    CHECK(split[0][0].cellid == 0);
    CHECK(split[0][1].cellid == 1);
    CHECK(split[1][0].cellid == 2);
    CHECK(split[3][4].cellid == 9);


    REQUIRE_THROWS(LoadBalancerBase::partition(problems, {3,3,3}));


}



TEST_CASE("LoadBalancerBase sendRecv() swap test"){

    
    using namespace Foam;

    
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

  

    DynamicList<DynamicList<ChemistryProblem>> send_buffer;
    send_buffer.setSize(1);
    send_buffer[0] = create_problems(10);

    auto recv_buffer = testableLoadBalancing::test_sendRecv<ChemistryProblem, Pstream::commsTypes::nonBlocking>(
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

TEST_CASE("LoadBalancerBase sendRecv() self send") {

    using namespace Foam;

    auto problems = create_problems(10);

    std::vector<int> sources      = {Pstream::myProcNo()};
    std::vector<int> destinations = {Pstream::myProcNo()};
    std::vector<int> counts       = {10};

    auto send_buffer = LoadBalancerBase::partition(problems, counts);

    auto recv_buffer = testableLoadBalancing::test_sendRecv<ChemistryProblem,
         Pstream::commsTypes::nonBlocking>(send_buffer, sources, destinations);


    CHECK(recv_buffer.size() == 1);
    
    auto recv_problems = recv_buffer[0];

    for (size_t i = 0; i < 10; ++i){
        CHECK(recv_problems[i].cellid == i);
    }


}

TEST_CASE("LoadBalancerBase sendRecv() self send 2") {

    using namespace Foam;

    auto problems = create_problems(1);

    std::vector<int> sources      = {Pstream::myProcNo()};
    std::vector<int> destinations = {Pstream::myProcNo()};
    std::vector<int> counts       = {1};

    auto send_buffer = LoadBalancerBase::partition(problems, counts);

    auto recv_buffer = testableLoadBalancing::test_sendRecv<ChemistryProblem,
         Pstream::commsTypes::nonBlocking>(send_buffer, sources, destinations);


    CHECK(recv_buffer.size() == 1);
    
    auto recv_problems = recv_buffer[0];

    for (size_t i = 0; i < 1; ++i){
        CHECK(recv_problems[i].cellid == i);
    }


}

/*



TEST_CASE("LoadBalancerBase balance() / unbalance()"){

    

    
    using namespace Foam;

    

    testableLoadBalancing l;

    auto problems = create_problems(50);
    l.updateState(problems);

    


    for (const auto& p : problems ){

    }

    //l.setState(l.determine_state(loads));

    label startOfRequests = Pstream::nRequests();
    auto buffer = l.balance(problems);

    if (Pstream::myProcNo() == Pstream::master()){
        CHECK(buffer.size() == Pstream::nProcs() - 1);
    }

    else{
        CHECK(buffer.size() == 0);
    }
        
    
    //Pstream::waitRequests(startOfRequests);

    if (buffer.size() > 0){
        auto buffer2 = l.unbalance(buffer);
    }

    
    
    
    
    



}

*/