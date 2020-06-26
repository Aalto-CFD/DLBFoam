#include "simpleBalancingMethod.H"

namespace Foam{


chemistryLoadBalancingMethod::sendRecvInfo simpleBalancingMethod::determine_state(const DynamicList<chemistryLoad>& loads) const{
    // This chunk balances the load between ranks n and n+1, periodically.

    // The ratio of (between 0, 1) original problems on a rank to swap between the neighbour.
    // zero means no problems are swapped and 1 means that all problems are swapped 
    float mutliplier = 0;

    size_t my_idx = rank_to_load_idx(loads, Pstream::myProcNo());
    auto my_load = loads[my_idx];

    
     
    
    if (Pstream::parRun()){


        if (Pstream::myProcNo() % 2 == 0){
        
            sendRecvInfo ret;
            

            int pair_rank = Pstream::myProcNo() + 1;
            if (pair_rank > Pstream::nProcs()){
                pair_rank = 0;
            } 
            auto load_pair = loads[rank_to_load_idx(loads, pair_rank)];
            int recv_count = int(mutliplier*load_pair.number_of_active_cells);
            int send_count = 0;
            int remaining_count = my_load.number_of_active_cells - send_count;
            
            ret.sources = {Pstream::myProcNo(), load_pair.rank }; 
            ret.destinations = {Pstream::myProcNo()};
            ret.number_of_problems = {remaining_count, recv_count};
            return ret;
        }


        else{
        
            sendRecvInfo ret;
            
            int pair_rank = Pstream::myProcNo() - 1;
            if (pair_rank < 0){
                pair_rank = Pstream::nProcs() - 1;
            }
            auto load_pair = loads[rank_to_load_idx(loads, pair_rank)];
            int recv_count = 0;
            int send_count = int(mutliplier*my_load.number_of_active_cells);
            int remaining_count = my_load.number_of_active_cells - send_count;
            
            ret.sources = {Pstream::myProcNo()};
            ret.destinations = {Pstream::myProcNo(),load_pair.rank }; 
            ret.number_of_problems = {remaining_count, send_count};
            return ret;
        }
    
    }
    

    //others do nothing
    sendRecvInfo ret;
    ret.sources = {Pstream::myProcNo()};
    ret.destinations = {Pstream::myProcNo()};
    ret.number_of_problems = {my_load.number_of_active_cells};
    

    return ret;
    
 


}

} //namespace Foam
