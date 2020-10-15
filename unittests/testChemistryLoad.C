#include "../third_party/catch.hpp"
#include "ChemistryLoad.H"
#include "ChemistryProblem.H"
#include "Pstream.H"




TEST_CASE("ChemistryLoad comparisons"){


    using namespace Foam;

    ChemistryLoad l1(1, 1.0);
    ChemistryLoad l2(1, 2.0);


    CHECK(!(l1 == l2));
    CHECK(l1 != l2);
    CHECK(l1 <= l2);
    CHECK(l2 >= l1);
    CHECK(l1 < l2);
    CHECK(l2 > l1);


}

TEST_CASE("ChemistryLoad mpi send and receive"){



    using namespace Foam;

    if (Pstream::parRun()) {



    ChemistryLoad l1(43, 1.433531);

    ChemistryLoad l2(-1, -2);


    int source = 0;
    int destination = 1;

    PstreamBuffers pBufs(Pstream::commsTypes::blocking);


    if (Pstream::myProcNo() == source){

        UOPstream send(destination, pBufs);
        send << l1;

    }

    pBufs.finishedSends();

    if (Pstream::myProcNo() == destination){
        UIPstream recv(source, pBufs);
        recv >> l2;
        CHECK(l1.rank == l2.rank);
        CHECK(l1.value == l2.value);
    }


    }



}