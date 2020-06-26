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

Base OF-dev git commit : 3c643bf2f8eb50b0531f40ba93be6ecb2364c955 
Base OF-dev file path : src/thermophysicalModels/chemistryModel/chemistrySolver/chemistrySolver/makePexChemistrySolverTypes.H

\*---------------------------------------------------------------------------*/

#ifndef makePexChemistrySolverTypes_H
#define makePexChemistrySolverTypes_H




#include "loadBalancedChemistryModel.H"
//#include "TDAC_loadBalancedChemistryModel.H" #TODO: add TDAC

#include "noChemistrySolver.H"
//#include "./ode/ode.H" //important to not use the stock ode for pyJac

#include "ode.H" //important to use the stock ODE's with others
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePexChemistrySolverType(SS, Comp, Thermo)                              \
                                                                               \
                                                                               \
    typedef SS<loadBalancedChemistryModel<Comp, Thermo>> loadBalanced##SS##Comp##Thermo;     \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        loadBalanced##SS##Comp##Thermo,                                               \
        (#SS"<" + word(loadBalancedChemistryModel<Comp, Thermo>::typeName_()) + "<"   \
        + word(Comp::typeName_()) + "," + Thermo::typeName() + ">>").c_str(),  \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<loadBalanced##SS##Comp##Thermo>               \
        add##loadBalanced##SS##Comp##Thermo##thermo##ConstructorTo##BasicChemistryModel\
##Comp##Table_;                                                                \
                                                                               \
                                                                               \


//TODO: add TDAC


#define makePexChemistrySolverTypes(Comp, Thermo)                                 \
                                                                               \
    makePexChemistrySolverType                                                    \
    (                                                                          \
        noChemistrySolver,                                                     \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makePexChemistrySolverType                                                    \
    (                                                                          \
        ode,                                                                   \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
