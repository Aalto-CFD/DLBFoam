/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Base OF-dev git commit : 4bd59d64db473a4f2a3654dc93081d958559dcad 
Base OF-dev file path : src/thermophysicalModels/chemistryModel/chemistrySolver/chemistrySolver/makePexChemistrySolvers.C

\*---------------------------------------------------------------------------*/

#include "makePexChemistrySolverTypes.H"

#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makePexChemistrySolverTypes(psiReactionThermo, constGasHThermoPhysics);
    makePexChemistrySolverTypes(psiReactionThermo, gasHThermoPhysics);

    makePexChemistrySolverTypes(rhoReactionThermo, constGasHThermoPhysics);
    makePexChemistrySolverTypes(rhoReactionThermo, gasHThermoPhysics);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
