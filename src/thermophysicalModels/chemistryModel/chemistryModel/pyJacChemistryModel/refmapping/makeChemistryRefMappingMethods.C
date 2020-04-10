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

\*---------------------------------------------------------------------------*/

#include "makeChemistryRefMappingMethods.H"

#include "thermoPhysicsTypes.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistryRefMappingMethods(psiReactionThermo, constGasHThermoPhysics);
    makeChemistryRefMappingMethods(psiReactionThermo, gasHThermoPhysics);


    makeChemistryRefMappingMethods
    (
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryRefMappingMethods
    (
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistryRefMappingMethods(psiReactionThermo, icoPoly8HThermoPhysics);
    makeChemistryRefMappingMethods(psiReactionThermo, constFluidHThermoPhysics);
    makeChemistryRefMappingMethods
    (
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistryRefMappingMethods(psiReactionThermo, constHThermoPhysics);

    makeChemistryRefMappingMethods(rhoReactionThermo, constGasHThermoPhysics);
    makeChemistryRefMappingMethods(rhoReactionThermo, gasHThermoPhysics);
    makeChemistryRefMappingMethods
    (
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryRefMappingMethods
    (
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistryRefMappingMethods(rhoReactionThermo, icoPoly8HThermoPhysics);
    makeChemistryRefMappingMethods(rhoReactionThermo, constFluidHThermoPhysics);
    makeChemistryRefMappingMethods
    (
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistryRefMappingMethods(rhoReactionThermo, constHThermoPhysics);

    // Chemistry solvers based on sensibleInternalEnergy

    makeChemistryRefMappingMethods(psiReactionThermo, constGasEThermoPhysics);
    makeChemistryRefMappingMethods(psiReactionThermo, gasEThermoPhysics);
    makeChemistryRefMappingMethods
    (
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryRefMappingMethods
    (
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistryRefMappingMethods(psiReactionThermo, icoPoly8EThermoPhysics);
    makeChemistryRefMappingMethods(psiReactionThermo, constFluidEThermoPhysics);
    makeChemistryRefMappingMethods
    (
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistryRefMappingMethods(psiReactionThermo, constEThermoPhysics);

    makeChemistryRefMappingMethods(rhoReactionThermo, constGasEThermoPhysics);
    makeChemistryRefMappingMethods(rhoReactionThermo, gasEThermoPhysics);
    makeChemistryRefMappingMethods
    (
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryRefMappingMethods
    (
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistryRefMappingMethods(rhoReactionThermo, icoPoly8EThermoPhysics);
    makeChemistryRefMappingMethods(rhoReactionThermo, constFluidEThermoPhysics);
    makeChemistryRefMappingMethods
    (
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistryRefMappingMethods(rhoReactionThermo, constEThermoPhysics);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
