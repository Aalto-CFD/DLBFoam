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
    struct sendRecvInfo {
        // would like to use DynamicList here but it appears it sucks ass and
        // nothing works as expected when its used
        // DynamicList<int> sources;             //ranks which send to this process
        // DynamicList<int> destinations;       //ranks to which this process sends to
        // DynamicList<int> number_of_problems; //number of problems which this rank sends/receivs

        std::vector<int> sources;            // ranks which send to this process
        std::vector<int> destinations;       // ranks to which this process sends to
        std::vector<int> number_of_problems; // number of problems which this rank sends/receivs
    };

public:
    // Constructors
    chemistryLoadBalancingMethod() = default;

    DynamicList<chemistryLoad> get_loads() const;

    DynamicList<DynamicList<chemistryProblem>>
    get_problems(const DynamicList<chemistryProblem>& all_problems) const {

        auto sources      = m_current_state.sources;
        auto destinations = m_current_state.destinations;
        auto counts       = m_current_state.number_of_problems;

        // TODO: Check only in debug mode
        if (sources.size() != 0 && destinations.size() != 0) {
            throw "A process can be either a source, destination or a do nothinger";
        }

        //TODO: replace with all to all
        auto send_buffers = get_send_buffer(all_problems);

        DynamicList<DynamicList<chemistryProblem>> recv_buffers;

        for (size_t i = 0; i < sources.size(); ++i){
            Isend_recv(send_buffers[i], recv_buffers[i], sources[i], counts[i]);
        }
        
        
        return recv_buffers;
    }

    DynamicList<DynamicList<chemistryProblem>>
    get_send_buffer(const DynamicList<chemistryProblem>& all_problems) const;

    DynamicList<chemistrySolution>
    get_solutions(const DynamicList<chemistrySolution>& current) const;

    ///
    ///@brief Similar to MPI_ISendrecv but for generic types with << and >> operators overloaded.
    ///       The communication is non-blocking.
    ///
    ///@tparam T Type with Ostream <<, Istream >> operators
    ///@param send_buffer the send data
    ///@param recv_buffer buffer to place the received data
    ///@param source source rank
    ///@param destination destination rank
    ///
    template <class T>
    void Isend_recv(const T& send_buffer, T& recv_buffer, int source, int destination) const {

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        if (Pstream::myProcNo() == source) {

            UOPstream send(destination, pBufs);
            send << send_buffer;
        }

        pBufs.finishedSends();

        if (Pstream::myProcNo() == destination) {
            UIPstream recv(source, pBufs);
            recv >> recv_buffer;
        }
    }

    void apply_balancing(const chemistryRefMappingMethod* mapper, PtrList<volScalarField>& Y) const;

    void print_parameters() const { std::cout << "I am balancing load" << std::endl; }

    bool active() const { return true; }

    //- Destructor
    virtual ~chemistryLoadBalancingMethod() {}

protected:
    void set_state(const sendRecvInfo& state) { m_current_state = state; }

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

private:
    sendRecvInfo m_current_state;
};

} // namespace Foam

#endif

// ************************************************************************* //
