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

#include <algorithm> //std::sort / std::partition_copy
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
    struct sendRecvInfo {

        std::vector<int> sources;            // ranks which send to this process
        std::vector<int> destinations;       // ranks to which this process sends to
        std::vector<int> number_of_problems; // number of problems which this rank sends/receivs
    };

    template <class T> using buffer_t = DynamicList<DynamicList<T>>;

private:
    sendRecvInfo m_current_state; // the current state of the object

public:
    void update_state(const DynamicList<chemistryProblem>& problems) {

        DynamicList<chemistryLoad> loads = get_loads(problems);
        set_state(determine_state(loads));
    }

    ///
    ///@brief Given a list of values splits them evenly between the mpi-processes.
    ///
    ///@tparam T type of value (typically a chemistryProblem)
    ///@param values values to split
    ///@return buffer_t<T> additional values coming to the current process
    ///
    template <class T> buffer_t<T> balance(const DynamicList<T>& values) const {

        // auto send_buffer = get_send_buffer(values);
        return send_recv<T, Pstream::commsTypes::blocking>(
            get_send_buffer(values), m_current_state.sources, m_current_state.destinations);

        // auto send_buffer = get_send_buffer(values);
        //return send_recv<T, Pstream::commsTypes::nonBlocking>(
        //    get_send_buffer(values), m_current_state.sources, m_current_state.destinations);    
    }

    ///
    ///@brief Given a buffer of values, sends the values back to their rightful owner processes.
    ///
    ///@tparam T type of value (typically a chemistrySolution)
    ///@param values values to send back
    ///@return buffer_t<T> values coming back to this process
    ///
    template <class T> buffer_t<T> unbalance(const buffer_t<T>& values) const {

        
        // It is probably impossible to avoid blocking in this phase
        return send_recv<T, Pstream::commsTypes::blocking>(
            values, m_current_state.destinations, m_current_state.sources);
    }


    
    template <class T> buffer_t<T> get_send_buffer(const DynamicList<T>& values) const {

        if (m_current_state.destinations.size() == 1) { return buffer_t<T>(1, values); }

        return partition(values, m_current_state.number_of_problems);
    }

    

    

    // Constructors
    chemistryLoadBalancingMethod() = default;

    ///
    ///@brief Check whether the load balancer is active or not.
    ///
    ///@return true
    ///@return false
    ///
    bool active() const { return true; }

    //- Destructor
    virtual ~chemistryLoadBalancingMethod() {}

    ///
    ///@brief Splits the big_array to arrays of size subsizes.
    ///
    ///@tparam T element type
    ///@param big_array the array to split
    ///@param subsizes  vector of subarray sizes
    ///@return DynamicList<DynamicList<T>> array of subarrays, if succesful
    ///
    template <class T>
    static DynamicList<DynamicList<T>> partition(const DynamicList<T>&   big_array,
                                                 const std::vector<int>& subsizes) {

        // TODO: only in debug mode
        if (std::accumulate(subsizes.begin(), subsizes.end(), 0) != big_array.size()) {

            int s_total = std::accumulate(subsizes.begin(), subsizes.end(), 0);
            int a_total = big_array.size();

            throw error(
                std::string(
                    "Total send count not matching the problem count for subsize total of: ") +
                std::to_string(s_total) + std::string(" and array total of: ") +
                std::to_string(a_total));
        }

        // allocate storage based on the first element
        T element_copy;
        if (big_array.size() > 0) { element_copy = big_array[0]; }

        DynamicList<DynamicList<T>> ret;

        for (auto s : subsizes) { ret.append(DynamicList<T>(s, element_copy)); }

        int start = 0;
        for (size_t i = 0; i < subsizes.size(); ++i) {
            int end = start + subsizes[i];
            std::copy(
                std::begin(big_array) + start, std::begin(big_array) + end, std::begin(ret[i]));
            start = end;
        }

        return ret;
    }

    ///
    ///@brief Sets the current state
    ///
    ///@param state state to set to
    ///
    void set_state(const sendRecvInfo& state) {

        m_current_state = state;

        if (!valid_state()) {

            print_state();

            throw error("Invalid balance state in chemistryLoadBalancing.");
        }
    }

    ///
    ///@brief Prints the current state information.
    ///
    ///
    void print_state() const {

        Pout << "Rank: " << Pstream::myProcNo()
             << " sends to: " << vector_to_string(m_current_state.destinations)
             << " receives from: " << vector_to_string(m_current_state.sources)
             << " counts: " << vector_to_string(m_current_state.number_of_problems) << endl;
    }

    template <class T> std::string vector_to_string(const std::vector<T>& vec) const {

        std::string ret = "(";

        for (const auto& v : vec) { ret += std::to_string(v) + ' '; }
        ret += ")";
        return ret;
    }

protected:
    size_t rank_to_load_idx(const DynamicList<chemistryLoad>& loads, int rank) const;

    bool valid_state() const {

        auto sources      = m_current_state.sources;
        auto destinations = m_current_state.destinations;

        if (sources.size() < 1) return false;
        if (destinations.size() < 1) return false;
        if (sources.size() > Pstream::nProcs()) return false;
        if (destinations.size() > Pstream::nProcs()) return false;
        if (sources[0] != Pstream::myProcNo())
            return false; // the first index must always be this process
        if (sources[0] != destinations[0]) return false;

        return true;
    }

    DynamicList<chemistryLoad> get_loads(const DynamicList<chemistryProblem>& problems) const;

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
    static buffer_t<ET> send_recv(const buffer_t<ET>&     send_buffer,
                                  const std::vector<int>& sources,
                                  const std::vector<int>& destinations) {

        if (Pstream::parRun()) {
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

        return send_buffer;
    }

private:
    ///
    ///@brief Given a list of problems belonging to this rank, determine the load of this process
    ///
    ///@param problems problems of this rank
    ///@return chemistryLoad my load
    ///
    virtual chemistryLoad get_my_load(const DynamicList<chemistryProblem>& problems) const;

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
