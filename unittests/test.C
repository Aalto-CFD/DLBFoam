#define CATCH_CONFIG_RUNNER

#include "../third_party/catch.hpp"

#include "IPstream.H"
#include "OPstream.H"
#include "UPstream.H"


int main(int argc, char* argv[])
{

    using namespace Foam;

//    if (Pstream::parRun()){
//        bool ok = UPstream::init(argc, argv, true);
//    }
    
    Catch::Session session;
  
    const int result = session.run(argc, argv);

    
    if (Pstream::parRun()){
        UPstream::exit(0);
    }
    //MPI_Finalize();

    return result;
    
    
}
