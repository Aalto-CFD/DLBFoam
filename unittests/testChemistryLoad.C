#include "../third_party/catch.hpp"
#include "chemistryLoad.H"
#include "Pstream.H"

TEST_CASE("chemistryLoad comparisons"){


    using namespace Foam;

    chemistryLoad l1(1, 1.0, 100);
    chemistryLoad l2(1, 2.0, 100);

    
    CHECK(!(l1 == l2));
    CHECK(l1 != l2);
    CHECK(l1 <= l2);
    CHECK(l2 >= l1);
    CHECK(l1 < l2);
    CHECK(l2 > l1);
        


}

TEST_CASE("chemistryLoad mpi send and receive"){


    using namespace Foam;

    chemistryLoad l1(43, 1.433531, 11330);

    chemistryLoad l2(-1, -2, -3);

    
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
        CHECK(l1.number_of_active_cells == l2.number_of_active_cells);
    }

    

    

    //fromNeighb >> nbrdAlpha >> nbrpFaceCells >> nbrp_rghPross >> nbrCellC ;


    


}