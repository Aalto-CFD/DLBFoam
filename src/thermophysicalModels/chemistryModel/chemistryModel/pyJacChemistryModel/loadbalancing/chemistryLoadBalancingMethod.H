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


#include <vector> //std::vector
#include <algorithm> //std::sort

#include "IOdictionary.H"
#include "scalarField.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"

#include "psiReactionThermo.H"

#include "chemistryProblem.H"
#include "chemistryLoad.H"
#include "chemistryRefMappingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//TODO: rename to something better
    enum WHATTODO {    e_SENDER = 0, 
                        e_RECEIVER = 1,
                        e_NEITHER = 2 
                    }; 
  



/*---------------------------------------------------------------------------*\
                           Class chemistryLoadBalancingMethod Declaration
\*---------------------------------------------------------------------------*/

class chemistryLoadBalancingMethod
{

public:

    
    // Constructors
    chemistryLoadBalancingMethod() = default;

    DynamicList<chemistryLoad> get_loads() const;


    void send_recv(const DynamicList<chemistryProblem>& problems, int source, int destination);

    





    void apply_balancing(const chemistryRefMappingMethod* mapper, PtrList<volScalarField>& Y) const;


    void print_parameters() const {
        std::cout << "I am balancing load" << std::endl;
    }

    bool active() const{return true;}

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
    virtual WHATTODO determine_state(const DynamicList<chemistryLoad>& loads) const = 0;


   
    




};

  


} //namespace Foam


#endif

// ************************************************************************* //
