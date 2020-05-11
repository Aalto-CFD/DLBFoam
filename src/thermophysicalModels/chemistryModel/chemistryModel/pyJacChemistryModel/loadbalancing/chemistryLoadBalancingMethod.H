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



    template<class T> using buffer_t = DynamicList<DynamicList<T>>;


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

    


    template<class T>
    buffer_t<T> balance(const DynamicList<T>& values) const{

    }


    ///
    ///@brief Given a list of problems splits them evenly between the mpi-processes.
    ///
    ///@param problems list of problems to split
    ///@return buffer_t<chemistryProblem> additional problems to solve
    ///
    buffer_t<chemistryProblem>
    get_problems(const DynamicList<chemistryProblem>& problems) const;


    ///
    ///@brief Returns the solutions back to the rightful owner processes.
    ///
    ///@param  solutions the solutions to send back to
    ///@return DynamicList<chemistrySolution> the solutions coming back to this process.
    ///
    DynamicList<chemistrySolution>
    get_solutions(const buffer_t<chemistrySolution>& solutions) const;


    


    ///
    ///@brief Sends the split send_buffer to sources and receives everything from destinations. Does not
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
    buffer_t<ET> send_recv(const buffer_t<ET>&      send_buffer,
                          const std::vector<int>& sources,
                          const std::vector<int>& destinations
                          ) const {

        

        PstreamBuffers pBufs(CT);

        for (int i = 0; i < destinations.size(); ++i){
            UOPstream send(destinations[i], pBufs);
            send << send_buffer[i];
            
        }
            
        pBufs.finishedSends();

        buffer_t<ET> recv_buffers;
        recv_buffers.setSize(sources.size());

        for (int i = 0; i < sources.size(); ++i){
            
            UIPstream recv(sources[i], pBufs);
            recv >> recv_buffers[i];
            
        }
        
        return recv_buffers;
        
    }


    ///
    ///@brief Given a list of problems splits the input list to sublists to be passed to other ranks
    ///
    ///@param all_problems list of problems that this rank has to solve
    ///@return DynamicList<DynamicList<chemistryProblem>> list of portions to be passed to other ranks
    ///
    DynamicList<DynamicList<chemistryProblem>>
    get_send_buffer(const DynamicList<chemistryProblem>& all_problems) const;

    

    
    //TODO: make private
    DynamicList<chemistryLoad> get_loads() const;
    




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
