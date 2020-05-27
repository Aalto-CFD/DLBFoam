#include "simpleLoadBalancing.H"

namespace Foam{


chemistryLoadBalancingMethod::sendRecvInfo simpleLoadBalancing::determine_state(const DynamicList<chemistryLoad>& loads) const{

    size_t my_idx = get_my_load_index(loads);
    auto my_load = loads[my_idx];

       
    if (my_idx == 0){
        //receive three problems from the most busy rank

        sendRecvInfo ret;
        int recv_count = 1;
        int send_count = 0;
        auto load_pair = loads[loads.size() - 1];
        int remaining_count = my_load.number_of_active_cells - send_count;
        
        ret.sources = {Pstream::myProcNo(), load_pair.rank }; //receive from busiest
        ret.destinations = {Pstream::myProcNo()};
        ret.number_of_problems = {remaining_count, recv_count};
        return ret;
    }

    else if (my_idx == loads.size() - 1){
        //send three problems to least busy
        sendRecvInfo ret;
        int recv_count = 0;
        int send_count = 1;
        auto load_pair = loads[0];
        int remaining_count = my_load.number_of_active_cells - send_count;
        
        ret.sources = {Pstream::myProcNo()};
        ret.destinations = {Pstream::myProcNo(),load_pair.rank }; //send to least busy
        ret.number_of_problems = {remaining_count, send_count};
        return ret;
    }
    



    //others do nothing
    sendRecvInfo ret;
    ret.sources = {Pstream::myProcNo()};
    ret.destinations = {Pstream::myProcNo()};
    ret.number_of_problems = {my_load.number_of_active_cells};

    Pout << Pstream::myProcNo() << endl;
    Pout << my_load.number_of_active_cells << endl;

    return ret;
    
 


}

} //namespace Foam
