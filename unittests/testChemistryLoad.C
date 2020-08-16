#include "../third_party/catch.hpp"
#include "chemistryLoad.H"
#include "chemistryProblem.H"
#include "Pstream.H"




TEST_CASE("chemistryLoad comparisons"){


    using namespace Foam;

    chemistryLoad l1(1, 1.0);
    chemistryLoad l2(1, 2.0);

    
    CHECK(!(l1 == l2));
    CHECK(l1 != l2);
    CHECK(l1 <= l2);
    CHECK(l2 >= l1);
    CHECK(l1 < l2);
    CHECK(l2 > l1);
    

    chemistryLoad l3(3, 1.0);

    l3+=1.5;
    l3-=1.0;
    l3*=2.0;
    l3/=2.0;

    CHECK(l3.rank == 3);
    CHECK(l3.value == 1.5);


    chemistryLoad l4(4, 3.0);
    chemistryLoad l5(5, 2.0);

    l4 /=l5;
    CHECK(l4.value == 1.5);
    CHECK(l4.value != 2.5);
    l4 += l5;
    l4 -=l5;
    l4 *=l5;
    CHECK(l4.rank == 4);
    CHECK(l4.value == 3.0);


}

TEST_CASE("chemistryLoad mpi send and receive"){

    

    using namespace Foam;

    if (Pstream::parRun()) {

    

    chemistryLoad l1(43, 1.433531);

    chemistryLoad l2(-1, -2);

    
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