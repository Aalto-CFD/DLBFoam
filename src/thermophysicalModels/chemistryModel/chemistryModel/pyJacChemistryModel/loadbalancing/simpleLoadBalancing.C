#include "simpleLoadBalancing.H"

namespace Foam{


bool simpleLoadBalancing::applyLoadBalancing() {
    //return check_if_refcell(Y,celli);
    return true;
}


List<scalarList> simpleLoadBalancing::compute_overhead(scalarField& t_cpu_list, labelField& ncells_list, labelField& nActiveCells_list)
{

    
    Info << "--Starting load balancing--" << endl;
    label nPs = Pstream::nProcs();
    label nVar = 5;

    
    // TODO: Check if these max, min, avg functions are accurate, if not revert to original implementation
    scalar maxCPUt = max(t_cpu_list);
    scalar minCPUt = min(t_cpu_list);
    scalar meanCPUt = average(t_cpu_list);  

    scalarList overheadT(nPs,0.0);
    scalarList nc_rm(nPs,0.0);  //- Relative cell count wrt. the average processor

    overheadT = t_cpu_list/meanCPUt;  //- TODO: Check if this makes sense, different from original

    Info << "Maximum, Minimum and Mean CPU use: "<< maxCPUt << ", " << minCPUt << ", " << meanCPUt << endl;


    //- If the mean CPU time is very low (e.g. too many cores), the minimal cpu times are
    //- large wrt. mean --> overhead vector is not realistic.
    //- Therefore we force the minimum element of the overhead vector to zero.
    forAll(overheadT, p_i)
    {
        if( (t_cpu_list[p_i] < 5*minCPUt) && (overheadT[p_i] < 0.9) )
        {
            overheadT[p_i] = 0.0;
        }
        // If the core has only a few chemically active cells, we don't want to include them into balancing
        if( (overheadT[p_i] >= 1.0) && (nActiveCells_list[p_i] < 50) )
        {
            overheadT[p_i] = 1.0; // not within tolerance --> not activated
        }    
        
        nc_rm[p_i] = ncells_list[p_i] * (1.0/max(overheadT[p_i],VSMALL));         
    }

    List<scalarList> overhead(2);
    overhead[0] = overheadT;
    overhead[1] = nc_rm;
    return overhead;
}

List<labelListList> simpleLoadBalancing::compute_global_stats(scalarField& t_cpu_list, labelField& nActiveCells_list, List<scalarList>& overhead)
{
    label nPs = Pstream::nProcs(); // TODO: Initialize these two in the constructor  
    label nVar = 5;
    scalarList overheadT = overhead[0];
    scalarList nc_rm = overhead[1];

    //- TODO: Make these arrays varying length not fixed
    labelList receiverInfo_mpi(mpiInfoTableSize_,0);
    labelListList receiverInfo_all(nPs);
    labelList senderInfo_mpi(mpiInfoTableSize_,0);
    labelListList senderInfo_all(nPs);


    for(int p_i=0; p_i<nPs; p_i++)
    {
        receiverInfo_all[p_i] = receiverInfo_mpi;
        senderInfo_all[p_i] = receiverInfo_mpi;
    } 

    int mpiBufferLim = mpiBufferLimit_; //- Mpi communication has an upper limit for data transfer. TODO: CHECK WITH HEIKKI
    scalar tol2bal = 0.02; //- Tolerance limit whether overhead is considered worth to balance or not
    scalar maxCPUt = max(t_cpu_list);

    //- If maxCPUt long enough, do balancing:
    if (maxCPUt > chemCPUTimeLimit_)
    {
        labelList order;
        sortedOrder(overheadT, order);
        sort(overheadT);

        for(int p_i = nPs-1; p_i>=0; p_i--) //- Iterate from largest to smallest to account for efficient filling of sent data
        {
            if(overheadT[p_i] >= (1.0+tol2bal))
            {
                int nCellsSentTot = 0;  //- Check that this does not go over nActiveCellsTot
                for(int s_i = 0; s_i<nPs; s_i++)  //- Iterate over all processors which we can send data
                {
                    scalar relFreeSpace = 1.0-overheadT[s_i];
                    scalar relData2Send = min(relFreeSpace,(overheadT[p_i]-1.0));  //- Amount of data between [0 1] to send
                    int nCells2Send = floor(relData2Send * nc_rm[order[p_i]]);

                    //- Prevent load balancing from sending more data than the core already has
                    int nCellsSentTot_tmp = nCellsSentTot + nCells2Send;
                    if(nCellsSentTot_tmp >= nActiveCells_list[order[p_i]])
                    {
                        nCells2Send = nCells2Send - (nCellsSentTot_tmp - nActiveCells_list[order[p_i]]) - 2; 
                    }
                    nCellsSentTot += nCells2Send;

                    //- If the following conditions apply, do load balancing

                    if((overheadT[p_i] > (1.0 + tol2bal)) && (relFreeSpace > tol2bal) && (s_i < p_i) && (nCells2Send>1))
                    {
                        overheadT[p_i] -= relData2Send;
                        overheadT[s_i] += relData2Send;

                        scalar mpiMessageSize = nVar*nCells2Send;
                        scalar mpiMessageSize_orig = mpiMessageSize;

                        //- Processor numbers and current list size for e.g. receiverInfo_all variable

                        int pid_r = order[s_i];
                        int pid_s = order[p_i];
                        int nLtmp_r = receiverInfo_all[pid_r][0];
                        int nLtmp_s = senderInfo_all[pid_s][0];

                        //- Following if/else defines whether you can send the data as one chunk or not
                        if(mpiMessageSize < mpiBufferLim)
                        {
                            receiverInfo_all[pid_r][nLtmp_r+1] = pid_s;    // processor id from which you receive information
                            senderInfo_all[pid_s][nLtmp_s+1] = pid_r;    // processor id to which you send information
                            int nctmp = nCells2Send;             
                            receiverInfo_all[pid_r][nLtmp_r+2] = nctmp;    
                            senderInfo_all[pid_s][nLtmp_s+2] = nctmp;
                            receiverInfo_all[pid_r][0] += 2;               // update table size
                            senderInfo_all[pid_s][0] += 2; 

                        }
                        else
                        {
                            bool dataAvailable = true; 
                            scalar mpiMessageSize_tmp = mpiMessageSize;
                            while (dataAvailable)
                            {
                                mpiMessageSize_tmp = min(mpiMessageSize,mpiBufferLim);
                                scalar nCells2Send_i = mpiMessageSize_tmp / mpiMessageSize_orig * nCells2Send;
                                int nctmp = floor(nCells2Send_i);              // number of cells to receive (note casting)
                                if( nctmp > 1 )    
                                {
                                    nLtmp_r = receiverInfo_all[pid_r][0];          // update list size because we are in while loop
                                    nLtmp_s = senderInfo_all[pid_s][0];      
                                
                                    receiverInfo_all[pid_r][nLtmp_r+1] = pid_s;    // processor id from which you receive information
                                    senderInfo_all[  pid_s][nLtmp_s+1] = pid_r;    // processor id to which you send information

                                    receiverInfo_all[pid_r][nLtmp_r+2] = nctmp;
                                    senderInfo_all[  pid_s][nLtmp_s+2] = nctmp;

                                    receiverInfo_all[pid_r][0] += 2;               // update table size
                                    senderInfo_all[pid_s][0] += 2;

                                }
                                mpiMessageSize -= mpiMessageSize_tmp; 

                                if(mpiMessageSize>0)
                                {
                                    dataAvailable = true; 
                                }    
                                else
                                {
                                    dataAvailable = false;
                                } 
                            }                            
                        }
                    }
                }
            }
        }
    }
    List<labelListList> global_stats(2);
    global_stats[0] = receiverInfo_all;
    global_stats[1] = senderInfo_all;
    return global_stats;

    //- Distibute the needed send/receive information to all processors
    /*
    for(int slave=Pstream::firstSlave(); slave<=Pstream::lastSlave(); slave++)        
    { 
        for (int j=0; j<mpiInfoTableSize_; j++)
        {
            receiverInfo[j] = receiverInfo_all[slave][j];
            senderInfo[j] = senderInfo_all[slave][j];
        } 
        OPstream toSlave(Pstream::commsTypes::blocking, slave);
        toSlave << receiverInfo;
        toSlave << senderInfo;
    }
    
    //Then the master itself
    int p_i_m = Pstream::masterNo();
    for (int j=0; j<mpiInfoTableSize_; j++)
    {
        receiverInfo[j] = receiverInfo_all[p_i_m][j];
        senderInfo[j] = senderInfo_all[p_i_m][j];
    }
    */    
}



} //namespace Foam