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

#include <vector>    //std::vector
#include <algorithm> //std::sort
#include <numeric>   //std::accumulate


#include "IOdictionary.H"
#include "scalarField.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"

#include "psiReactionThermo.H"

#include "chemistryProblem.H"
#include "chemistrySolution.H"
#include "chemistryLoad.H"
#include "chemistryRefMappingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                           Class chemistryLoadBalancingMethod Declaration
\*---------------------------------------------------------------------------*/

class chemistryLoadBalancingMethod
{

public:
    struct sendRecvInfo
    {
        DynamicList<int> sources;             //ranks which send to this process
        DynamicList<int> destinations;       //ranks to which this process sends to
        DynamicList<int> number_of_problems; //number of problems which this rank sends/receivs
    };

public:
    // Constructors
    chemistryLoadBalancingMethod() = default;

    DynamicList<chemistryLoad> get_loads() const;

    DynamicList<DynamicList<chemistryProblem>> get_problems(const DynamicList<chemistryProblem> &current) const{

        auto sources = m_current_state.sources;
        auto destinations = m_current_state.destinations;

        //TODO: Check only in debug mode
        if (sources.size() != 0 && destinations.size() != 0){
            throw "A process can be either a source, destination or a do nothinger";
        }

        auto send_buffer = get_send_buffer(current);

        /*
        DynamicList<DynamicList<chemistryProblem>> received_problems;

        for (const auto& source : sources){
            DynamicList<Problem> recv_buffer;


        }
        */
        return send_buffer;
        //return received_problems;


    }


    DynamicList<DynamicList<chemistryProblem>> get_send_buffer(const DynamicList<chemistryProblem>& current) const{

        int n_problems     = current.size();  
        
        int total_send_count =  std::accumulate(m_current_state.number_of_problems.begin(), 
                                                m_current_state.number_of_problems.end(),
                                                0);

        //number of own problems this process will solve
        int n_my_problems = n_problems - total_send_count;

        DynamicList<DynamicList<chemistryProblem>> send_buffer;

        int start = n_my_problems;
        for (const auto& n : m_current_state.number_of_problems){

            
            int end   = start + n;
            DynamicList<chemistryProblem> slice;
            std::copy(current.begin() + start, current.begin() + end, slice.begin());
            //DynamicList<chemistryProblem> slice(current.begin() + start, current.begin() + end);
            send_buffer.append(slice);
            start     = end; 
            

        }
        return send_buffer;

    }


    DynamicList<chemistrySolution> get_solutions(const DynamicList<chemistrySolution> &current) const;

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
    void Isend_recv(const T &send_buffer, T &recv_buffer, int source, int destination) const
    {

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        if (Pstream::myProcNo() == source)
        {

            UOPstream send(destination, pBufs);
            send << send_buffer;
        }

        pBufs.finishedSends();

        if (Pstream::myProcNo() == destination)
        {
            UIPstream recv(source, pBufs);
            recv >> recv_buffer;
        }
    }

    void apply_balancing(const chemistryRefMappingMethod *mapper, PtrList<volScalarField> &Y) const;

    void print_parameters() const
    {
        std::cout << "I am balancing load" << std::endl;
    }

    bool active() const { return true; }

    //- Destructor
    virtual ~chemistryLoadBalancingMethod() {}

protected:
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
    virtual sendRecvInfo determine_state(const DynamicList<chemistryLoad> &loads) const = 0;




private:
 
    sendRecvInfo m_current_state;



};

} //namespace Foam

#endif

// ************************************************************************* //
