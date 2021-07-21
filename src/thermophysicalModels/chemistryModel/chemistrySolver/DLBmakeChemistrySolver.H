/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing 
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           | 
     \\/     M anipulation  | 2020, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of DLBFoam library, derived from OpenFOAM.

    https://github.com/Aalto-CFD/DLBFoam

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


#define defineChemistrySolver(Model, ThermoPhysics)                            \
                                                                               \
    typedef Model<ThermoPhysics>                                               \
        Model##ThermoPhysics;                                                  \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Model##ThermoPhysics,                                                  \
        (                                                                      \
            word(Model##ThermoPhysics::typeName_())                            \
          + "<" + ThermoPhysics::typeName() + ">"                              \
        ).c_str(),                                                             \
        0                                                                      \
    )


#define makeChemistrySolver(Solver, Model, ThermoPhysics)                      \
                                                                               \
    typedef Solver<Model<ThermoPhysics>>                                       \
        Solver##Model##ThermoPhysics;                                          \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Solver##Model##ThermoPhysics,                                          \
        (                                                                      \
            word(Solver##Model##ThermoPhysics::typeName_())                    \
          + "<" + word(Model<ThermoPhysics>::typeName_())                      \
          + "<" + ThermoPhysics::typeName() + ">>"                             \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        basicChemistryModel,                                                   \
        Solver##Model##ThermoPhysics,                                          \
        thermo                                                                 \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistrySolvers(nullArg, ThermoPhysics)                         \
    defineChemistrySolver                                                      \
    (                                                                          \
        LoadBalancedChemistryModel,                                            \
        ThermoPhysics                                                          \
    );                                                                         \
                                                                               \
    defineChemistrySolver                                                      \
    (                                                                          \
        pyJacLoadBalancedChemistryModel,                                       \
        ThermoPhysics                                                          \
    )


#define makeChemistrySolvers(Solver, ThermoPhysics)                            \
    makeChemistrySolver                                                        \
    (                                                                          \
        Solver,                                                                \
        LoadBalancedChemistryModel,                                            \
        ThermoPhysics                                                          \
    );                                                                         \
    makeChemistrySolver                                                        \
    (                                                                          \
        Solver,                                                                \
        pyJacLoadBalancedChemistryModel,                                       \
        ThermoPhysics                                                          \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
