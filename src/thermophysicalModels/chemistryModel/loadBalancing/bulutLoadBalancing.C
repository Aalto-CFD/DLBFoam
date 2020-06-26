#include "bulutLoadBalancing.H"


namespace Foam{


    

chemistryLoadBalancingMethod::sendRecvInfo 
bulutLoadBalancing::determine_state(const DynamicList<chemistryLoad>& loads) const {
    //- Assumes load array in the shape loads[i] = chemCPUT[i]/meanCPUT
    //- Assumes that corrections (e.g. too idle ranks) are already made        
    // TODO: Implement a minCPUTime limit to start considering load balancing        
    sendRecvInfo curr_state;
    
    //this is unnecessary, please formulate this function so that no copy of the loads is needed
    //TODO: reformulate the function to slightly more sane version
    DynamicList<chemistryLoad> loads_copy = loads; 


    // iterate from largest to smallest to account for efficient filling of the sent data
    for(int sender_idx=Pstream::nProcs()-1; sender_idx>=0; sender_idx--) 
    {

        chemistryLoad& sender_load = loads_copy[sender_idx];


        if(sender_load.value >= (1.0+m_tolerance_to_balance))
        {
            int nCellsSentTot = 0;
            // iterate over all ranks which we can send data
            for(int receiver_idx=0; receiver_idx<sender_idx; receiver_idx++) 
            {

                chemistryLoad& receiver_load = loads_copy[receiver_idx];


                int nCells2Send = ncells_to_send(sender_load, receiver_load, nCellsSentTot);
                nCellsSentTot +=  nCells2Send;                         
                if( large_sender(sender_load)  && small_receiver(receiver_load) )
                {    
                    update_load( sender_load, receiver_load );

                    if(sender_load.rank  == Pstream::myProcNo())
                    {
                        curr_state.destinations.push_back(receiver_load.rank);
                        curr_state.number_of_problems.push_back(nCells2Send);
                    }
                    else if(receiver_load.rank == Pstream::myProcNo())
                    {
                        curr_state.sources.push_back(sender_load.rank);
                        curr_state.number_of_problems.push_back(nCells2Send);
                    }
                }           
            }
        }
    }


    
    size_t my_idx = rank_to_load_idx(loads, Pstream::myProcNo());
    auto my_load = loads[my_idx];

    
    int send_problems = std::accumulate(curr_state.number_of_problems.begin(), curr_state.number_of_problems.end(), 0);
    int remaining_problems = my_load.number_of_active_cells - send_problems;


    curr_state.destinations.insert(curr_state.destinations.begin(), Pstream::myProcNo());
    curr_state.sources.insert(curr_state.sources.begin(), Pstream::myProcNo());
    curr_state.number_of_problems.insert(curr_state.number_of_problems.begin(), remaining_problems);
    return curr_state;
}



void bulutLoadBalancing::update_load(chemistryLoad& sender_load, chemistryLoad& receiver_load) const
{
    sender_load.value -= min(1.0 - receiver_load.value, sender_load.value - 1.0);
    receiver_load.value += min(1.0 - receiver_load.value, sender_load.value - 1.0);
}



int bulutLoadBalancing::ncells_to_send(chemistryLoad& sender_load, chemistryLoad& receiver_load, int nCellsSentTot) const
{
    scalar relFreeSpace = 1.0 - receiver_load.value;
    scalar relData2send = min(relFreeSpace,(sender_load.value-1.0));  
    scalar nc_rm = sender_load.number_of_active_cells * (1.0/max(sender_load.value,1e-12));
    int nCells2Send =  floor(relData2send * nc_rm);

    int nCellsSentTot_tmp = nCellsSentTot + nCells2Send;
    if( nCellsSentTot_tmp >= sender_load.number_of_active_cells) 
    {
        nCells2Send = nCells2Send - (nCellsSentTot_tmp - sender_load.number_of_active_cells) - 2;
    }
    return nCells2Send; 
}





bool bulutLoadBalancing::large_sender(const chemistryLoad& load) const
{
    if(load.value>1.0+m_tolerance_to_balance)
    {
        return true;
    }
    return false;
}

bool bulutLoadBalancing::small_receiver(const chemistryLoad& load) const
{
    if(load.value < 1.0-m_tolerance_to_balance)
    {
        return true;
    }
    return false;
}







}