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
#include "chemistrySolution.H"

#ifdef DEBUG
static inline void runtime_assert(bool condition, const char* msg) {
    if (!condition) throw std::runtime_error(msg);
}
#else
static inline void runtime_assert([[maybe_unused]] bool condition, [[maybe_unused]] const char* msg) {}
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                           Class chemistryLoadBalancingMethod Declaration
\*---------------------------------------------------------------------------*/

class chemistryLoadBalancingMethod {

public:
    struct sendRecvInfo {

        std::vector<label> sources;            // ranks which send to this process
        std::vector<label> destinations;       // ranks to which this process sends to
        std::vector<label> number_of_problems; // number of problems which this rank sends/receivs
    };

    template <class T> using buffer_t = DynamicList<DynamicList<T>>;

private:
    sendRecvInfo m_current_state; // the current state of the object

public:


    chemistryLoadBalancingMethod() = default;

    virtual ~chemistryLoadBalancingMethod() {}


    static double get_mean(const DynamicList<chemistryLoad>& loads);
    static chemistryLoad get_min(const DynamicList<chemistryLoad>& vec);
    static chemistryLoad get_max(const DynamicList<chemistryLoad>& vec);

    chemistryLoad
    compute_my_load(const DynamicList<chemistryProblem>& problems) {
        auto   lambda = [](scalar sum, const chemistryProblem& rhs) { return sum + rhs.cpuTime; };
        scalar sum    = std::accumulate(problems.begin(), problems.end(), scalar(0), lambda);
        return chemistryLoad(Pstream::myProcNo(), sum);
    }

    template<class T>
    static DynamicList<T> all_gather(const T& my_data) {


        DynamicList<T> ret(Pstream::nProcs(), my_data);
        ret[Pstream::myProcNo()] = my_data;

        label tag                    = 1;
        Pstream::gatherList(ret, tag);
        Pstream::scatterList(ret, tag);


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
    ///@brief Get the state object
    ///
    ///@return const sendRecvInfo& 
    ///
    const sendRecvInfo& get_state() const {
        return m_current_state;
    }



    virtual void update_state(const DynamicList<chemistryProblem>& problems) = 0;

    /*
    void update_state(const DynamicList<chemistryProblem>& problems) {

        DynamicList<chemistryLoad> loads = get_loads(problems);
        set_state(determine_state(loads));
    }
    */

    ///
    ///@brief Given a list of values splits them evenly between the mpi-processes.
    ///
    ///@tparam T type of value (typically a chemistryProblem)
    ///@param values values to split
    ///@return buffer_t<T> additional values coming to the current process
    ///
    template <class T> buffer_t<T> balance(const DynamicList<T>& values) const {

    
        //return send_recv<T, Pstream::commsTypes::blocking>(
        //    get_send_buffer(values), m_current_state.sources, m_current_state.destinations);

        
        return send_recv<T, Pstream::commsTypes::nonBlocking>(
            get_send_buffer(values), m_current_state.sources, m_current_state.destinations);    
    }

    ///
    ///@brief Given a buffer of values, sends the values back to their rightful owner processes.
    ///
    ///@tparam T type of value (typically a chemistrySolution)
    ///@param values values to send back
    ///@return buffer_t<T> values coming back to this process
    ///
    template <class T> buffer_t<T> unbalance(const buffer_t<T>& values) const {

        
        //return send_recv<T, Pstream::commsTypes::blocking>(
        //    values, m_current_state.destinations, m_current_state.sources);

        // It is probably impossible to avoid blocking in this phase
        return send_recv<T, Pstream::commsTypes::nonBlocking>(
            values, m_current_state.destinations, m_current_state.sources);

    }


    
    template <class T> buffer_t<T> get_send_buffer(const DynamicList<T>& values) const {

        if (m_current_state.destinations.size() == 1) { return buffer_t<T>(1, values); }

        return partition(values, m_current_state.number_of_problems);
    }

    

    


    ///
    ///@brief Check whether the load balancer is active or not.
    ///
    ///@return true
    ///@return false
    ///
    bool active() const { return true; }


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
        
        runtime_assert(check_sendcount_size(big_array,subsizes), "Total send count not matching the problem count!");

        // allocate storage based on the first element
        T element_copy;
        if (big_array.size() > 0) { element_copy = big_array[0]; }

        DynamicList<DynamicList<T>> ret;

        for (auto s : subsizes) { ret.append(DynamicList<T>(s, element_copy)); }

        label start = 0;
        for (label i = 0; i < subsizes.size(); ++i) {
            label end = start + subsizes[i];
            std::copy(
                std::begin(big_array) + start, std::begin(big_array) + end, std::begin(ret[i]));
            start = end;
        }

        return ret;
    }

    ///
    ///@brief Checks if the total send count matches the problem count for subsizes
    ///
    ///@tparam T element type
    ///@param big_array the array to split
    ///@param subsizes  vector of subarray sizes
    ///@return bool
    ///
    template <class T>
    static bool check_sendcount_size(const DynamicList<T>&   big_array,
                              const std::vector<label>& subsizes) {
        if (std::accumulate(subsizes.begin(), subsizes.end(), 0) != big_array.size()) {
            /*
            int s_total = std::accumulate(subsizes.begin(), subsizes.end(), 0);
            int a_total = big_array.size();

            throw error(
                std::string(
                    "Total send count not matching the problem count for subsize total of: ") +
                std::to_string(s_total) + std::string(" and array total of: ") +
                std::to_string(a_total));
            */
            return false;
        } 
            return true;
    }

    

    ///
    ///@brief Prints the current state information.
    ///
    ///
    void print_state() const {

        //receiver
        if (m_current_state.sources.size() > 1){
            Pout << "Receiver rank: " << Pstream::myProcNo()
                 << " receives from "  << vector_to_string(m_current_state.sources)
                 << " own problems "   << vector_to_string(m_current_state.number_of_problems)
                 << endl;

        }
        //sender
        else if(m_current_state.destinations.size() > 1) {
            Pout << "Sender rank: " << Pstream::myProcNo()
             << " sends to: " << vector_to_string(m_current_state.destinations)
             << " counts: " << vector_to_string(m_current_state.number_of_problems) 
             << " remaining problems:  " << m_current_state.number_of_problems.back() << endl;
        }

        else {
            Pout << "Rank: " << Pstream::myProcNo() << " does not take part in balancing.";
        }
        
    }

    template <class T> std::string vector_to_string(const std::vector<T>& vec) const {

        std::string ret = "(";

        for (const auto& v : vec) { ret += std::to_string(v) + ' '; }
        ret += ")";
        return ret;
    }

protected:
    label rank_to_load_idx(const DynamicList<chemistryLoad>& loads, label rank) const;

    bool valid_state() const {

        /*
        auto sources      = m_current_state.sources;
        auto destinations = m_current_state.destinations;

        if (sources.size() < 1) return false;
        if (destinations.size() < 1) return false;
        if (sources.size() > size_t(Pstream::nProcs())) return false;
        if (destinations.size() > size_t(Pstream::nProcs())) return false;
        if (sources[0] != Pstream::myProcNo())
            return false; // the first index must always be this process
        if (sources[0] != destinations[0]) return false;
        */
        //TODO: enable
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
                                  const std::vector<label>& sources,
                                  const std::vector<label>& destinations) {

        if (Pstream::parRun()) {
            PstreamBuffers pBufs(CT);

            for (label i = 0; i < destinations.size(); ++i) {
                UOPstream send(destinations[i], pBufs);
                send << send_buffer[i];
            }
            pBufs.finishedSends();

            buffer_t<ET> recv_buffers;
            recv_buffers.setSize(sources.size());
            for (label i = 0; i < sources.size(); ++i) {
                UIPstream recv(sources[i], pBufs);
                recv >> recv_buffers[i];
            }

            return recv_buffers;
        }

        return send_buffer;
    }




};

} // namespace Foam

#endif

// ************************************************************************* //
