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

#include "./ode_pyJac/ode_pyJac.H"
#include "LoadBalancedChemistryModel.H"
#include "pyJacLoadBalancedChemistryModel.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

#include "forCommonGases.H"
#include "forCommonLiquids.H"
#include "forPolynomials.H"
#include "DLBmakeChemistrySolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCommonGases(makeChemistrySolvers, ode_pyJac, psiReactionThermo);
    forCommonGases(makeChemistrySolvers, ode_pyJac, rhoReactionThermo);

    forCommonLiquids(makeChemistrySolvers, ode_pyJac, rhoReactionThermo);

    forPolynomials(makeChemistrySolvers, ode_pyJac, rhoReactionThermo);
}


// ************************************************************************* //
