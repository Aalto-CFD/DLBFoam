/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::chemistryLoadBalancingMethod

Description
    An abstract class for chemistry tabulation.

SourceFiles
    chemistryLoadBalancingMethod.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryLoadBalancingMethod_H
#define chemistryLoadBalancingMethod_H

#include <algorithm> //std::sort
#include <numeric>   //std::accumulate
#include <vector>    //std::vector

#include "IOdictionary.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

#include "psiReactionThermo.H"

#include "chemistryLoad.H"
#include "chemistryProblem.H"
#include "chemistryRefMappingMethod.H"
#include "chemistrySolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                           Class chemistryLoadBalancingMethod Declaration
\*---------------------------------------------------------------------------*/

class chemistryLoadBalancingMethod {



public:
    template <class T> using buffer_t = DynamicList<DynamicList<T>>;

    struct sendRecvInfo {

        std::vector<int> sources;            // ranks which send to this process
        std::vector<int> destinations;       // ranks to which this process sends to
        std::vector<int> number_of_problems; // number of problems which this rank sends/receivs
    };

private:
    sendRecvInfo m_current_state; // the current state of the object

    //TODO: make user input
    static constexpr scalar m_tolerance_to_balance = 0.02;


public:



    // Constructors
    chemistryLoadBalancingMethod() = default;



    void update_state() {
        //- Assumes load array in the shape loads[i] = chemCPUT[i]/meanCPUT
        //- Assumes that corrections (e.g. too idle ranks) are already made        
        // TODO: Implement a minCPUTime limit to start considering load balancing        
        sendRecvInfo curr_state;
        DynamicList<chemistryLoad> loads  = get_loads();
        // iterate from largest to smallest to account for efficient filling of the sent data
        for(int sender_idx=Pstream::nProcs()-1; sender_idx>=0; sender_idx--) 
        {
            if(loads[sender_idx].value >= (1.0+m_tolerance_to_balance))
            {
                int nCellsSentTot = 0;
                // iterate over all ranks which we can send data
                for(int receiver_idx=0; receiver_idx<sender_idx; receiver_idx++) 
                {
                    int nCells2Send = ncells_to_send(loads, sender_idx, receiver_idx, nCellsSentTot);
                    nCellsSentTot +=  nCells2Send;                         
                    if( large_sender(loads[sender_idx])  && small_receiver(loads[receiver_idx]) )
                    {    
                        update_load( loads[sender_idx], loads[receiver_idx] );

                        if(loads[sender_idx].rank  == Pstream::myProcNo())
                        {
                            curr_state.destinations.push_back(loads[receiver_idx].rank);
                            curr_state.number_of_problems.push_back(nCells2Send);
                        }
                        else if(loads[receiver_idx].rank == Pstream::myProcNo())
                        {
                            curr_state.sources.push_back(loads[sender_idx].rank);
                            curr_state.number_of_problems.push_back(nCells2Send);
                        }
                    }           
                }
            }
        }        
        set_state(curr_state);
    }

    

    void update_load(chemistryLoad& sender_load, chemistryLoad& receiver_load) const
    {
        sender_load.value -= min(1.0 - receiver_load.value, sender_load.value - 1.0);
        receiver_load.value += min(1.0 - receiver_load.value, sender_load.value - 1.0);
    }

    int ncells_to_send(const DynamicList<chemistryLoad>& loads, const int& sender_idx, const int& receiver_idx, const int& nCellsSentTot)
    {
        scalar relFreeSpace = 1.0 - loads[receiver_idx].value;
        scalar relData2send = min(relFreeSpace,(loads[sender_idx].value-1.0));  
        scalar nc_rm = loads[sender_idx].number_of_active_cells * (1.0/max(loads[sender_idx].value,1e-12));
        int nCells2Send =  floor(relData2send * nc_rm);

        int nCellsSentTot_tmp = nCellsSentTot + nCells2Send;
        if( nCellsSentTot_tmp >= loads[sender_idx].number_of_active_cells) 
        {
            nCells2Send = nCells2Send - (nCellsSentTot_tmp - loads[sender_idx].number_of_active_cells) - 2;
        }
        return nCells2Send; 
    }

    bool large_sender(const chemistryLoad& load) const
    {
        if(load.value>1.0+m_tolerance_to_balance)
        {
            return true;
        }
        return false;
    }

    bool small_receiver(const chemistryLoad& load) const
    {
        if(load.value < 1.0-m_tolerance_to_balance)
        {
            return true;
        }
        return false;
    }


    ///
    ///@brief Given a list of values splits them evenly between the mpi-processes.
    ///
    ///@tparam T type of value (typically a chemistryProblem)
    ///@param values values to split
    ///@return buffer_t<T> additional values coming to the current process
    ///
    template <class T> buffer_t<T> balance(const DynamicList<T>& values) const {

        auto sources      = m_current_state.sources;
        auto destinations = m_current_state.destinations;

        // TODO: Check only in debug mode
        if (sources.size() != 0 && destinations.size() != 0) {
            throw "A process can be either a source, destination or a do nothinger";
        }

        auto send_buffer = get_send_buffer(values);

        // The unittests currently test such that everyone sends to master which requires a blocking
        // commstype in reality we dont want to block at this phase
        // TODO: change to nonBlocking
        return send_recv<T, Pstream::commsTypes::blocking>(send_buffer, sources, destinations);
    }

    ///
    ///@brief Given a buffer of values, sends the values back to their rightful owner processes.
    ///
    ///@tparam T type of value (typically a chemistrySolution)
    ///@param values values to send back
    ///@return buffer_t<T> values coming back to this process
    ///
    template <class T> buffer_t<T> unbalance(const buffer_t<T>& values) const {

        auto sources      = m_current_state.sources;
        auto destinations = m_current_state.destinations;

        // TODO: Check only in debug mode
        if (sources.size() != 0 && destinations.size() != 0) {
            throw "A process can be either a source, destination or a do nothinger";
        }

        // It is probably impossible to avoid blocking in this phase
        return send_recv<T, Pstream::commsTypes::blocking>(values, destinations, sources);
    }



    ///
    ///@brief Check whether the load balancer is active or not.
    ///
    ///@return true 
    ///@return false 
    ///
    bool active() const { return true; }

    //- Destructor
    virtual ~chemistryLoadBalancingMethod() {}





protected:
    ///
    ///@brief Sets the current state
    ///
    ///@param state state to set to
    ///
    void set_state(const sendRecvInfo& state) { m_current_state = state; }

    ///
    ///@brief Gathers the load information from all processes
    ///
    ///@return DynamicList<chemistryLoad> list size nprocesses of load objects
    ///
    DynamicList<chemistryLoad> get_loads() const;

    ///
    ///@brief Given a list of problems splits the input list to sublists to be passed to other ranks
    ///
    ///@param all_problems list of problems that this rank has to solve
    ///@return DynamicList<DynamicList<chemistryProblem>> list of portions to be passed to other
    /// ranks
    ///
    DynamicList<DynamicList<chemistryProblem>>
    get_send_buffer(const DynamicList<chemistryProblem>& all_problems) const;

    ///
    ///@brief Sends the split send_buffer to sources and receives everything from destinations. Does
    /// not
    ///       allow for overlapping sources or destinations.
    ///
    ///@tparam ET element type of the send and receive buffers
    ///@tparam CT mpi communication type blocking/nonBlocking
    ///@param send_buffer the buffer to send
    ///@param sources the ranks to receive from
    ///@param destinations the ranks to send to
    ///@return buffer_t<ET> the received buffer
    ///
    template <class ET, Pstream::commsTypes CT>
    buffer_t<ET> send_recv(const buffer_t<ET>&     send_buffer,
                           const std::vector<int>& sources,
                           const std::vector<int>& destinations) const {

        PstreamBuffers pBufs(CT);

        for (int i = 0; i < destinations.size(); ++i) {
            UOPstream send(destinations[i], pBufs);
            send << send_buffer[i];
        }
        pBufs.finishedSends();

        buffer_t<ET> recv_buffers;
        recv_buffers.setSize(sources.size());
        for (int i = 0; i < sources.size(); ++i) {
            UIPstream recv(sources[i], pBufs);
            recv >> recv_buffers[i];
        }

        return recv_buffers;
    }

private:
    ///
    ///@brief Determines the load of this process
    ///
    ///@return chemistryLoad load of this process
    ///
    virtual chemistryLoad get_my_load() const = 0;

    ///
    ///@brief Determines if a process is a sender receiver or do-nothinger
    ///
    ///@param loads a sorted array of loads
    ///@return WHATTODO state information on what the process should do
    ///
    virtual sendRecvInfo determine_state(const DynamicList<chemistryLoad>& loads) const = 0;
};

} // namespace Foam

#endif

// ************************************************************************* //
