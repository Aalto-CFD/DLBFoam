/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           |
     \\/     M anipulation  | 2022, Aalto University, Finland
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



#if ${method}CppTest == chemistryModelCppTest

#include "makeChemistryReductionMethod.H"

namespace Foam
{
    defineChemistryReductionMethod(nullArg, ThermoPhysics);
}

#include "noChemistryReduction.H"
#include "DAC.H"
#include "DRG.H"
#include "DRGEP.H"
#include "EFA.H"
#include "PFA.H"

namespace Foam
{
    makeChemistryReductionMethod(none, ThermoPhysics);
    makeChemistryReductionMethod(DAC, ThermoPhysics);
    makeChemistryReductionMethod(DRG, ThermoPhysics);
    makeChemistryReductionMethod(DRGEP, ThermoPhysics);
    makeChemistryReductionMethod(EFA, ThermoPhysics);
    makeChemistryReductionMethod(PFA, ThermoPhysics);
}

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "makeReaction.H"

#include "ArrheniusReactionRate.H"
#include "LandauTellerReactionRate.H"
#include "thirdBodyArrheniusReactionRate.H"

#include "JanevReactionRate.H"
#include "powerSeriesReactionRate.H"

#include "LangmuirHinshelwoodReactionRate.H"
#include "MichaelisMentenReactionRate.H"

#include "ChemicallyActivatedReactionRate.H"
#include "FallOffReactionRate.H"

#include "LindemannFallOffFunction.H"
#include "SRIFallOffFunction.H"
#include "TroeFallOffFunction.H"

namespace Foam
{
    defineReaction(nullArg, ThermoPhysics);

    makeIRNReactions(ArrheniusReactionRate, ThermoPhysics);
    makeIRNReactions(LandauTellerReactionRate, ThermoPhysics);
    makeIRNReactions(thirdBodyArrheniusReactionRate, ThermoPhysics);
    makeIRReactions(JanevReactionRate, ThermoPhysics);
    makeIRReactions(powerSeriesReactionRate, ThermoPhysics);

    makeIRReactions(LangmuirHinshelwoodReactionRate, ThermoPhysics);
    makeIReactions(MichaelisMentenReactionRate, ThermoPhysics);

    makeIRRPressureDependentReactions
    (
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction,
        ThermoPhysics
    );

    makeIRRPressureDependentReactions
    (
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction,
        ThermoPhysics
    );

    makeIRRPressureDependentReactions
    (
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction,
        ThermoPhysics
    );

    makeIRRPressureDependentReactions
    (
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction,
        ThermoPhysics
    );

    makeIRRPressureDependentReactions
    (
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction,
        ThermoPhysics
    );

    makeIRRPressureDependentReactions
    (
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction,
        ThermoPhysics
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluxLimitedLangmuirHinshelwoodReactionRate.H"

namespace Foam
{
    makeGeneralReaction
    (
        IrreversibleReaction,
        fluxLimitedLangmuirHinshelwoodReactionRate,
        ThermoPhysics
    );
}


// ************************************************************************* //
