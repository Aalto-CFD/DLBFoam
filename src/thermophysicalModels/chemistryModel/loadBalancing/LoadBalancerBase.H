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
    Foam::LoadBalancerBase

Description
    An abstract base class for all balancing methods.

SourceFiles
    LoadBalancerBase.C

\*---------------------------------------------------------------------------*/

#ifndef LoadBalancerBase_H
#define LoadBalancerBase_H

#include "ChemistryLoad.H"
#include "ChemistryProblem.H"
#include "ChemistrySolution.H"
#include "runtime_assert.H"

#include <algorithm> //std::min/max element
#include <numeric>   //std::accumulate
#include <vector>    //std::vector

namespace Foam
{

class LoadBalancerBase
{

public:


    struct BalancerState{
        std::vector<label> sources; //ranks which send to this process
        std::vector<label> destinations; //ranks to which this process sends to
        std::vector<label> nProblems; //counts which this process sends/receivs
        label nRemaining; //own problems / solutions

    };


    template <class T>
    using buffer_t = DynamicList<DynamicList<T>>;

private:
    BalancerState state_; // the current state of the object

public:
    LoadBalancerBase() = default;

    virtual ~LoadBalancerBase() = default;

public:
    ///
    ///@brief The load balancing algorithm which each derived class must define
    ///
    ///@param problems Problems of this process
    ///
    virtual void updateState(const DynamicList<ChemistryProblem>& problems) = 0;



    ///
    ///@brief Computes the mean of the input loads
    ///
    ///@param loads Loads to compute the means from
    ///@return double mean load
    ///
    static double getMean(const DynamicList<ChemistryLoad>& loads);

    ///
    ///@brief Finds the minimum of the input loads
    ///
    ///@param loads Loads to find the min from
    ///@return ChemistryLoad Minimum load
    ///
    static ChemistryLoad getMin(const DynamicList<ChemistryLoad>& loads);

    ///
    ///@brief Finds the maximum of the input loads
    ///
    ///@param loads Loads to find the max from
    ///@return ChemistryLoad Maximum load
    ///
    static ChemistryLoad getMax(const DynamicList<ChemistryLoad>& loads);

    ///
    ///@brief Computes the load based on the input problems
    ///
    ///@param problems Problems to compute the load from
    ///@return ChemistryLoad The (rank specific) load of the list of problems
    ///
    static ChemistryLoad
    computeLoad(const DynamicList<ChemistryProblem>& problems);

    ///
    ///@brief Gathers the data from all processes
    ///
    ///@tparam T a serializable type
    ///@param myData Data from this process
    ///@return DynamicList<T> List of datas from all processes with List[i]
    /// corresponding to the data from rank i
    ///
    template <class T>
    static DynamicList<T> allGather(const T& myData);

    ///
    ///@brief Sets the current state to the input state
    ///
    ///@param state State to set to
    ///
    void setState(const BalancerState& state);

    ///
    ///@brief Get the current state
    ///
    ///@return const BalancerState& current state
    ///
    const BalancerState& getState() const
    {
        return state_;
    }




    ///
    ///@brief Given a list of values splits them evenly between the
    /// mpi-processes.
    ///
    ///@tparam T type of value (typically a ChemistryProblem)
    ///@param values values to split
    ///@return buffer_t<T> additional values coming to the current process
    ///
    template <class T>
    buffer_t<T> balance(const DynamicList<T>& values) const;

    ///
    ///@brief Given a buffer of values, sends the values back to their rightful
    /// owner processes.
    ///
    ///@tparam T type of value (typically a ChemistrySolution)
    ///@param values values to send back
    ///@return buffer_t<T> values coming back to this process
    ///
    template <class T>
    buffer_t<T> unbalance(const buffer_t<T>& values) const;

    ///
    ///@brief Get the send buffer based on the current state and the input
    ///values
    ///
    ///@tparam T Any serializable type
    ///@param values Values to put to the send buffer
    ///@return buffer_t<T> The send buffer
    ///
    template <class T>
    buffer_t<T> getSendBuffer(const DynamicList<T>& values) const;

    ///
    ///@brief Check whether the load balancer is active or not.
    ///
    ///@return true
    ///@return false
    ///
    bool active() const;

    ///
    ///@brief Splits the bigArray to arrays of size subsizes.
    ///
    ///@tparam T element type
    ///@param bigArray the array to split
    ///@param subsizes  vector of subarray sizes
    ///@return DynamicList<DynamicList<T>> array of subarrays, if succesful
    ///
    template <class T>
    static DynamicList<DynamicList<T>>
    partition(const DynamicList<T>& bigArray, const std::vector<int>& subsizes);

    ///
    ///@brief Checks if the total send count matches the problem count for
    /// subsizes
    ///
    ///@tparam T element type
    ///@param bigArray the array to split
    ///@param subsizes  vector of subarray sizes
    ///@return bool
    ///
    template <class T>
    static bool checkSizes(
        const DynamicList<T>& bigArray, const std::vector<label>& subsizes);

    ///
    ///@brief Prints the current state information.
    ///
    ///
    void printState() const;

    ///
    ///@brief Converts a vector to std::string for printing
    ///
    ///@tparam type satisfying std::to_string(T)
    ///@param vec Vector to convert
    ///@return std::string String representation of the vector
    ///
    template <class T>
    static std::string vectorToString(const std::vector<T>& vec);

    ///
    ///@brief Finds the index of the load having a load.rank = rank
    ///
    ///@param loads Loads to find the rank from
    ///@param rank  Rank to find
    ///@return label Index of the rank
    ///
    label
    rankToLoadIdx(const DynamicList<ChemistryLoad>& loads, label rank) const;

    ///
    ///@brief
    ///
    ///@return true
    ///@return false
    ///
    bool validState() const;

    ///
    ///@brief Sends the split send_buffer to sources and receives everything
    /// from destinations.
    ///
    ///@tparam ET element type of the send and receive buffers
    ///@tparam CT mpi communication type blocking/nonBlocking
    ///@param send_buffer the buffer to send
    ///@param sources the ranks to receive from
    ///@param destinations the ranks to send to
    ///@return buffer_t<ET> the received buffer
    ///
    template <class ET, Pstream::commsTypes CT>
    static buffer_t<ET> sendRecv(
        const buffer_t<ET>&       send_buffer,
        const std::vector<label>& sources,
        const std::vector<label>& destinations);



    
    template <class ET>
    static void sendCall(
        const buffer_t<ET>&       send_buffer,
        const std::vector<label>& destinations,
        PstreamBuffers& pBufs)
    {

        if(Pstream::parRun())
        {
            for(label i = 0; i < label(destinations.size()); ++i)
            {
                UOPstream send(destinations[i], pBufs);
                send << send_buffer[i];
            }
            //pBufs.finishedSends();
        }    


    }


    template <class ET>
    static void recvCall(
        buffer_t<ET>&       recv_buffer,
        const std::vector<label>& sources,
        PstreamBuffers& pBufs)
    {
        if(Pstream::parRun()){

            //pBufs.finishedSends();
             
            recv_buffer.setSize(sources.size());
            for(label i = 0; i < label(sources.size()); ++i)
            {
                UIPstream recv(sources[i], pBufs);
                recv >> recv_buffer[i];
            }
        }

    }




    template<class T>
    DynamicList<T> getRemaining(const DynamicList<T>& values) {
    
        DynamicList<T> remaining(state_.nRemaining, T());
    
        std::copy(values.end() - state_.nRemaining,
             values.end(),
             remaining.begin()
        );

        return remaining;
    }


};

template <class T>
DynamicList<T> LoadBalancerBase::allGather(const T& myData)
{

    DynamicList<T> ret(Pstream::nProcs(), myData);
    ret[Pstream::myProcNo()] = myData;

    label tag = 1;
    Pstream::gatherList(ret, tag);
    Pstream::scatterList(ret, tag);

    return ret;
}

template <class T>
LoadBalancerBase::buffer_t<T>
LoadBalancerBase::balance(const DynamicList<T>& values) const
{

    auto sendBuffer = partition(values, state_.nProblems);

    auto recvBuffer = sendRecv<T, Pstream::commsTypes::nonBlocking>(
        sendBuffer, state_.sources, state_.destinations);

    return recvBuffer;
    

    //return sendRecv<T, Pstream::commsTypes::nonBlocking>(
    //    getSendBuffer(values), state_.sources, state_.destinations);
}





template <class T>
LoadBalancerBase::buffer_t<T>
LoadBalancerBase::unbalance(const buffer_t<T>& values) const
{
    
     
    return sendRecv<T, Pstream::commsTypes::nonBlocking>(
        values, state_.destinations, state_.sources);
}

template <class T>
LoadBalancerBase::buffer_t<T>
LoadBalancerBase::getSendBuffer(const DynamicList<T>& values) const
{
    return partition(values, state_.nProblems);
}

template <class T>
LoadBalancerBase::buffer_t<T>
LoadBalancerBase::partition(
    const DynamicList<T>& bigArray, const std::vector<int>& subsizes)
{

    //TODO: dont make copies, use subList instead

    buffer_t<T> ret;

    label start = 0;
    for (const auto& count : subsizes){
        label end = start + count;
        DynamicList<T> slice(count, T());
        std::copy(
            std::begin(bigArray) + start,
            std::begin(bigArray) + end,
            std::begin(slice)
        );
        ret.append(slice);
        start = end;
    }


    return ret;


}

template <class T>
bool LoadBalancerBase::checkSizes(
    const DynamicList<T>& bigArray, const std::vector<label>& subsizes)
{
    if(std::accumulate(subsizes.begin(), subsizes.end(), 0) != label(bigArray.size()))
    {
        return false;
    }
    return true;
}

template <class T>
std::string
LoadBalancerBase::vectorToString(const std::vector<T>& vec)
{

    std::string ret = "(";

    for(const auto& v : vec)
    {
        ret += std::to_string(v) + ' ';
    }
    ret += ")";
    return ret;
}

template <class ET, Pstream::commsTypes CT>
LoadBalancerBase::buffer_t<ET>
LoadBalancerBase::sendRecv(
    const LoadBalancerBase::buffer_t<ET>& send_buffer,
    const std::vector<label>&                         sources,
    const std::vector<label>&                         destinations)
{

    if(Pstream::parRun())
    {

    
        PstreamBuffers pBufs(CT);

        for(label i = 0; i < label(destinations.size()); ++i)
        {
            UOPstream send(destinations[i], pBufs);
            send << send_buffer[i];
        }
        pBufs.finishedSends();

        buffer_t<ET> recv_buffers;
        recv_buffers.setSize(sources.size());
        for(label i = 0; i < label(sources.size()); ++i)
        {
            UIPstream recv(sources[i], pBufs);
            recv >> recv_buffers[i];
        }

        return recv_buffers;
    
    }

    return send_buffer;
}

} // namespace Foam

#endif

// ************************************************************************* //
