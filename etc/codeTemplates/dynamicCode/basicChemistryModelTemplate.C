/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           |
     \\/     M anipulation  | 2022-2023, Aalto University, Finland
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

#include "typedefThermo.H"

#include "${specie}.H"

#include "thermo.H"

// EoS
#include "${equationOfState}.H"

// Thermo
#include "${thermo}Thermo.H"
#include "${energy}.H"

// Transport
#include "${transport}Transport.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ${SHA1sum}
    //
    // Unique function name that can be checked if the correct library version
    // has been loaded
    void ${typeName}_${SHA1sum}(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define ThermoPhysics                                                          \
    ${transport}Transport${energy}${thermo}Thermo${equationOfState}${specie}


namespace Foam
{
    typedefThermo
    (
        ${transport}Transport,
        ${energy},
        ${thermo}Thermo,
        ${equationOfState},
        ${specie}
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define chemistryModelCppTest 0
#define loadBalancedChemistryModelCppTest 1
#define loadBalanced_pyJacChemistryModelCppTest 2

#include "makeChemistrySolver.H"
#if ${method}ChemistryModelCppTest == loadBalancedChemistryModelCppTest
    #include "${method}ChemistryModel.H"
#elif ${method}ChemistryModelCppTest == loadBalanced_pyJacChemistryModelCppTest
    #include "${method}ChemistryModel.H"
#else
    #include "${method}.H"
#endif
#include "${solver}.H"

namespace Foam
{
    #if ${method}ChemistryModelCppTest == loadBalanced_pyJacChemistryModelCppTest
    defineChemistrySolver(chemistryModel, ThermoPhysics);
    makeChemistrySolver(${solver}, chemistryModel, ThermoPhysics);
    defineChemistrySolver(loadBalancedChemistryModel, ThermoPhysics);
    makeChemistrySolver(${solver}, loadBalancedChemistryModel, ThermoPhysics);
    defineChemistrySolver(loadBalanced_pyJacChemistryModel, ThermoPhysics);
    makeChemistrySolver(${solver}, loadBalanced_pyJacChemistryModel, ThermoPhysics);
    #elif ${method}ChemistryModelCppTest == loadBalancedChemistryModelCppTest
    defineChemistrySolver(chemistryModel, ThermoPhysics);
    makeChemistrySolver(${solver}, chemistryModel, ThermoPhysics);
    defineChemistrySolver(loadBalancedChemistryModel, ThermoPhysics);
    makeChemistrySolver(${solver}, loadBalancedChemistryModel, ThermoPhysics);
    #elif ${method}CppTest == chemistryModelCppTest
    defineChemistrySolver(chemistryModel, ThermoPhysics);
    makeChemistrySolver(${solver}, chemistryModel, ThermoPhysics);
    #else
    defineChemistrySolver(${method}, ThermoPhysics);
    makeChemistrySolver(${solver}, ${method}, ThermoPhysics);
    #endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define chemistryModelMethod 0

#if ${method}Method == chemistryModelMethod

#include "makeChemistryReductionMethod.H"

#include "noChemistryReduction.H"
#include "DAC.H"
#include "DRG.H"
#include "DRGEP.H"
#include "EFA.H"
#include "PFA.H"

namespace Foam
{
    defineChemistryReductionMethod(nullArg, ThermoPhysics);
    
    makeChemistryReductionMethod(none, ThermoPhysics);
    makeChemistryReductionMethod(DAC, ThermoPhysics);
    makeChemistryReductionMethod(DRG, ThermoPhysics);
    makeChemistryReductionMethod(DRGEP, ThermoPhysics);
    makeChemistryReductionMethod(EFA, ThermoPhysics);
    makeChemistryReductionMethod(PFA, ThermoPhysics);
}

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define chemistryModelMethod 0

#if ${method}Method == chemistryModelMethod

#include "makeReaction.H"

#include "ArrheniusReactionRate.H"
#include "LandauTellerReactionRate.H"
#include "thirdBodyArrheniusReactionRate.H"

#include "JanevReactionRate.H"
#include "powerSeriesReactionRate.H"

#include "FallOffReactionRate.H"
#include "ChemicallyActivatedReactionRate.H"
#include "LindemannFallOffFunction.H"
#include "SRIFallOffFunction.H"
#include "TroeFallOffFunction.H"

#include "MichaelisMentenReactionRate.H"
#include "LangmuirHinshelwoodReactionRate.H"
#include "fluxLimitedLangmuirHinshelwoodReactionRate.H"
#include "surfaceArrheniusReactionRate.H"

namespace Foam
{
    defineReaction(nullArg, ThermoPhysics);

    // Irreversible/reversible/non-equilibrium-reversible reactions
    makeIRNReactions(ArrheniusReactionRate, ThermoPhysics);
    makeIRNReactions(LandauTellerReactionRate, ThermoPhysics);
    makeIRNReactions(thirdBodyArrheniusReactionRate, ThermoPhysics);

    // Irreversible/reversible reactions
    makeIRReactions(JanevReactionRate, ThermoPhysics);
    makeIRReactions(powerSeriesReactionRate, ThermoPhysics);

    // Irreversible/reversible fall-off reactions
    makeIRTemplate2Reactions
    (
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction,
        ThermoPhysics
    );
    makeIRTemplate2Reactions
    (
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction,
        ThermoPhysics
    );
    makeIRTemplate2Reactions
    (
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction,
        ThermoPhysics
    );

    // Irreversible/reversible chemically activated reactions
    makeIRTemplate2Reactions
    (
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction,
        ThermoPhysics
    );
    makeIRTemplate2Reactions
    (
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction,
        ThermoPhysics
    );
    makeIRTemplate2Reactions
    (
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction,
        ThermoPhysics
    );

    // Michaelis-Menten Reactions
    makeIReactions(MichaelisMentenReactionRate, ThermoPhysics);

    // Langmuir-Hinshelwood Reactions
    makeIRReactions(LangmuirHinshelwoodReactionRate, ThermoPhysics);

    // Flux-limited Langmuir-Hinshelwood Reactions
    makeGeneralReaction
    (
        IrreversibleReaction,
        fluxLimitedLangmuirHinshelwoodReactionRate,
        ThermoPhysics
    );

    // Surface-Arrhenius Reactions
    makeGeneralReaction
    (
        IrreversibleReaction,
        surfaceArrheniusReactionRate,
        ThermoPhysics
    );
}

#endif


// ************************************************************************* //
