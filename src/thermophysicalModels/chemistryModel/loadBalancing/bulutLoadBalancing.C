#include "bulutLoadBalancing.H"


namespace Foam{


    

void bulutLoadBalancing::update_state(const DynamicList<chemistryProblem>& problems) {


    auto my_load = compute_my_load(problems);
    auto loads = all_gather(my_load);

    std::sort(loads.begin(), loads.end()); 

    //sendRecvInfo ret = initialize_to_nothinger(my_load); 
    sendRecvInfo ret;
    // Only balance if maxCPUT is larger than user-defined value
    // TODO: make chemCPUTimeLimit user defined
    auto   sum_op = [](scalar sum, const chemistryLoad& rhs) { return sum + rhs.value; };
    scalar meanCPUT =
        std::accumulate(loads.begin(), loads.end(), scalar(0), sum_op) / loads.size();
    Info<<"MEANCPUT IS: "<<meanCPUT<<endl;
    scalar chemCPUTimeLimit = 2;
    label send_problems = 0;
    if(active())
    {   
        // iterate from largest to smallest to account for efficient filling of the sent data
        for(label sender_idx=Pstream::nProcs()-1; sender_idx>=0; sender_idx--) 
        {
            chemistryLoad& sender_load = loads[sender_idx];
            label index = 0;

            if(sender_load.value > (meanCPUT+m_tolerance_to_balance))
            {
                label nCellsSentTot = 0;
                // iterate over all ranks which we can send data
                for(label receiver_idx=0; receiver_idx<sender_idx; receiver_idx++) 
                {
                    chemistryLoad& receiver_load = loads[receiver_idx];
                    if( large_sender(sender_load,meanCPUT)  && small_receiver(receiver_load,meanCPUT) )
                    {
                        scalar load_to_send = cpu_time_to_send(sender_load, receiver_load, meanCPUT);
                        update_load(sender_load,receiver_load,load_to_send);
                        label nCells2Send;
                        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);                    

                        if(sender_load.rank  == Pstream::myProcNo())
                        {
                            nCells2Send = cpu_time_to_ncells(load_to_send,problems,index);
                            ret.destinations.push_back(receiver_load.rank);
                            ret.number_of_problems.push_back(nCells2Send);
                            nCellsSentTot +=  nCells2Send;
                            
                            UOPstream send(receiver_load.rank, pBufs);
                            send << nCells2Send;  
                        }
                                          
                        pBufs.finishedSends();  

                        if(receiver_load.rank == Pstream::myProcNo())
                        {

                            UIPstream recv(sender_load.rank, pBufs);
                            recv >> nCells2Send;
                            ret.sources.push_back(sender_load.rank);
                            ret.number_of_problems.push_back(nCells2Send);
                        }
                    }

                }
            }
        }
        send_problems = std::accumulate(ret.number_of_problems.begin(), ret.number_of_problems.end(), 0);

    }
    ret.sources.push_back(Pstream::myProcNo());
    ret.destinations.push_back(Pstream::myProcNo());
    if(ret.destinations.size()>1)
    {
        ret.number_of_problems.push_back(problems.size()-send_problems);
    }
    else
    {
        ret.number_of_problems.push_back(problems.size());
    }

    set_state(ret);


}

chemistryLoadBalancingMethod::sendRecvInfo 
bulutLoadBalancing::initialize_to_nothinger(const DynamicList<chemistryProblem>& problems) const
{
    sendRecvInfo ret;
    ret.sources = {Pstream::myProcNo()};
    ret.destinations = {Pstream::myProcNo()};
    ret.number_of_problems = {problems.size()};
    return ret;
}

void bulutLoadBalancing::update_self_send(sendRecvInfo& ret) const
{  
    if(ret.destinations.size()>1)
    {
        label send_problems = std::accumulate(ret.number_of_problems.begin()+1, ret.number_of_problems.end(), 0);
        ret.number_of_problems[0] -= send_problems;
    }
}



scalar bulutLoadBalancing::get_mean_load(DynamicList<chemistryLoad>& loads) const
{
    scalar meanCPUT = 0;
    forAll(loads,i)
    {
        meanCPUT+=loads[i].value/loads.size();
    }
    meanCPUT = max(meanCPUT,VSMALL);
    
    return meanCPUT;

}



scalar bulutLoadBalancing::cpu_time_to_send(chemistryLoad& sender_load, chemistryLoad& receiver_load, scalar& meanCPUT) const
{
    scalar loadToSend = min((meanCPUT-receiver_load.value),(sender_load.value-meanCPUT));
    return loadToSend;
}

void bulutLoadBalancing::update_load(chemistryLoad& sender_load, chemistryLoad& receiver_load, scalar& loadToSend) const
{
    sender_load.value -= loadToSend;
    receiver_load.value += loadToSend;
}


label bulutLoadBalancing::cpu_time_to_ncells(const scalar& load_to_send, const DynamicList<chemistryProblem>& all_problems, label& index) const
{
    label nCells2Send = 0;
    scalar relSendLoad = 0;
    for (label i  = index; i<all_problems.size();i++)
    {
        relSendLoad += all_problems[i].cpuTime;
        nCells2Send++;
        index++;
        if (relSendLoad>0.99*load_to_send)
        {
            break;
        }
    } 
    return nCells2Send;
}




bool bulutLoadBalancing::large_sender(const chemistryLoad& load, const scalar& meanCPUT) const
{
    if(load.value>meanCPUT)
    {
        return true;
    }
    return false;
}

bool bulutLoadBalancing::small_receiver(const chemistryLoad& load, const scalar& meanCPUT) const
{
    if(load.value<meanCPUT)
    {
        return true;
    }
    return false;
}







}