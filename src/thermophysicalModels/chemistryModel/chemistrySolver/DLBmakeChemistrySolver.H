/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing 
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           | 
     \\/     M anipulation  | 2020, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of DLBFoam library, derived from OpenFOAM.

    https://github.com/blttkgl/DLBFoam

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
    
\*---------------------------------------------------------------------------*/

#ifndef DLBmakeChemistrySolver_H
#define DLBmakeChemistrySolver_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistrySolver(Model, ReactionThermo, ThermoPhysics)            \
                                                                               \
    typedef Model<ReactionThermo, ThermoPhysics>                               \
        Model##ReactionThermo##ThermoPhysics;                                  \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Model##ReactionThermo##ThermoPhysics,                                  \
        (                                                                      \
            word(Model##ReactionThermo##ThermoPhysics::typeName_()) + "<"      \
          + ReactionThermo::typeName + "," + ThermoPhysics::typeName() + ">"   \
        ).c_str(),                                                             \
        0                                                                      \
    );


#define makeChemistrySolver(Solver, Model, ReactionThermo, ThermoPhysics)      \
                                                                               \
    typedef Solver<Model<ReactionThermo, ThermoPhysics>>                       \
        Solver##Model##ReactionThermo##ThermoPhysics;                          \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Solver##Model##ReactionThermo##ThermoPhysics,                          \
        (                                                                      \
            word(Solver##Model##ReactionThermo##ThermoPhysics::typeName_())    \
          + "<" + word(Model<ReactionThermo, ThermoPhysics>::typeName_())      \
          + "<" + word(ReactionThermo::typeName_()) + ","                      \
          + ThermoPhysics::typeName() + ">>"                                   \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicChemistryModel<ReactionThermo>::                                      \
        addthermoConstructorToTable                                            \
        <Solver##Model##ReactionThermo##ThermoPhysics>                         \
        add##Solver##Model##ReactionThermo##ThermoPhysics##\
thermoConstructorToTable_;

#define defineChemistrySolvers(ReactionThermo, ThermoPhysics)                  \
    defineChemistrySolver                                                      \
    (                                                                          \
        LoadBalancedChemistryModel,                                            \
        ReactionThermo,                                                        \
        ThermoPhysics                                                          \
    );                                                                         \
    defineChemistrySolver                                                      \
    (                                                                          \
        pyJacLoadBalancedChemistryModel,                                       \
        ReactionThermo,                                                        \
        ThermoPhysics                                                          \
    )

#define makeChemistrySolvers(Solver, ReactionThermo, ThermoPhysics)            \
    makeChemistrySolver                                                        \
    (                                                                          \
        Solver,                                                                \
        LoadBalancedChemistryModel,                                            \
        ReactionThermo,                                                        \
        ThermoPhysics                                                          \
    );                                                                         \
    makeChemistrySolver                                                        \
    (                                                                          \
        Solver,                                                                \
        pyJacLoadBalancedChemistryModel,                                       \
        ReactionThermo,                                                        \
        ThermoPhysics                                                          \
    )                                                               

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
