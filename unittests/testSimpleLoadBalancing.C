#include "../third_party/catch.hpp"

#include "simpleLoadBalancing.H"
#include "chemistryProblem.H"


namespace Foam{


//create some arbitrary data
DynamicList<chemistryProblem> create_problems2(int count){



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



TEST_CASE("simpleLoadBalancing balance() / unbalance()"){

    
    simpleLoadBalancing l;

    auto problems = create_problems2(3 + 1*Pstream::myProcNo());
    l.update_state(problems);


    auto buffer_s = l.balance(problems);

    Pout << buffer_s.size() << endl;

    //CHECK(buffer_s.size() != 0);


}



}