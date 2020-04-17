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


#include "IOdictionary.H"
#include "scalarField.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"

#include "psiReactionThermo.H"

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


    // Constructors
    chemistryLoadBalancingMethod() = default;

    


    virtual void apply_balancing(const chemistryRefMappingMethod* mapper, PtrList<volScalarField>& Y) const = 0;


    virtual bool active()           const = 0;
    virtual void print_parameters() const = 0;

    //- Destructor
    virtual ~chemistryLoadBalancingMethod() {}



    
   
    




};

  


} //namespace Foam


#endif

// ************************************************************************* //
