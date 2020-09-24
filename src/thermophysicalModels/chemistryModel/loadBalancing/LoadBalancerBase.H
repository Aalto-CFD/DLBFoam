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
#include "RecvBuffer.H"
#include "SendBuffer.H"
#include "runtime_assert.H"

#include <algorithm> //std::min/max element
#include <numeric>   //std::accumulate
#include <vector>    //std::vector

namespace Foam
{

class LoadBalancerBase
{

public:
    struct BalancerState
    {
        std::vector<label> sources;      // ranks which send to this process
        std::vector<label> destinations; // ranks to which this process sends to
        std::vector<label> nProblems; // counts which this process sends/receivs
        label nRemaining;             // own problems / solutions
    };

private:
    BalancerState state_; // the current state of the object

public:
    LoadBalancerBase() = default;

    virtual ~LoadBalancerBase() = default;

public:
    ///
    ///@brief Check whether load balancing is active or not
    ///
    ///@return true
    ///@return false
    ///
    virtual bool active() const = 0;

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
    ///@return RecvBuffer<T> additional values coming to the current process
    ///
    template <class T>
    RecvBuffer<T> balance(const DynamicList<T>& values) const;

    ///
    ///@brief Given a buffer of values, sends the values back to their rightful
    /// owner processes.
    ///
    ///@tparam T type of value (typically a ChemistrySolution)
    ///@param values values to send back
    ///@return RecvBuffer<T> values coming back to this process
    ///
    template <class T>
    RecvBuffer<T> unbalance(const RecvBuffer<T>& values) const;

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
    ///@tparam Indexable an indexable send buffer type
    ///@param send_buffer the buffer to send
    ///@param sources the ranks to receive from
    ///@param destinations the ranks to send to
    ///@return RecvBuffer<ET> the received buffer
    ///
    template <class ET, class Indexable>
    static RecvBuffer<ET> sendRecv(
        const Indexable&          send_buffer,
        const std::vector<label>& sources,
        const std::vector<label>& destinations);

    ///
    ///@brief Slices an nRemaining size portion from the _end_ of the values
    ///
    ///@tparam T element type of the input list
    ///@param values the input list to slice
    ///@return SubList<T> the sliced list
    ///
    template <class T>
    SubList<T> getRemaining(const DynamicList<T>& values)
    {

        return SubList<T>(
            values, state_.nRemaining, values.size() - state_.nRemaining);
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
RecvBuffer<T> LoadBalancerBase::balance(const DynamicList<T>& values) const
{

    return sendRecv<T, SendBuffer<T>>(
        SendBuffer<T>(values, state_.nProblems),
        state_.sources,
        state_.destinations);
}

template <class T>
RecvBuffer<T> LoadBalancerBase::unbalance(const RecvBuffer<T>& values) const
{

    return sendRecv<T, RecvBuffer<T>>(
        values, state_.destinations, state_.sources);
}

template <class T>
std::string LoadBalancerBase::vectorToString(const std::vector<T>& vec)
{

    std::string ret = "(";

    for(const auto& v : vec)
    {
        ret += std::to_string(v) + ' ';
    }
    ret += ")";
    return ret;
}

template <class ET, class Indexable>
RecvBuffer<ET> LoadBalancerBase::sendRecv(
    const Indexable&          send_buffer,
    const std::vector<label>& sources,
    const std::vector<label>& destinations)
{

    RecvBuffer<ET> ret;

    if(Pstream::parRun())
    {

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        for(label i = 0; i < label(destinations.size()); ++i)
        {
            UOPstream send(destinations[i], pBufs);
            send << send_buffer[i];
        }

        pBufs.finishedSends();

        ret.setSize(sources.size());
        for(label i = 0; i < label(sources.size()); ++i)
        {
            UIPstream recv(sources[i], pBufs);
            recv >> ret[i];
        }
    }

    return ret;
}

} // namespace Foam

#endif

// ************************************************************************* //
