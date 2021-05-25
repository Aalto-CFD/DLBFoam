#include "catch.hpp"

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


} //namespace Foam





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


    using send_buffer_t = DynamicList<DynamicList<ChemistryProblem>>;

    send_buffer_t send_buffer;
    send_buffer.setSize(1);
    send_buffer[0] = create_problems(10);

    auto recv_buffer = LoadBalancerBase::sendRecv<ChemistryProblem, send_buffer_t>(
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

