#define CATCH_CONFIG_RUNNER

#include "../third_party/catch.hpp"



int main(int argc, char* argv[])
{

    
    
    Catch::Session session;
  
    const int result = session.run(argc, argv);

    //Communication::Mpi::finalize();
    return result;
    
    
}
