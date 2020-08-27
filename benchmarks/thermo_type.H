#pragma once

//https://github.com/OpenFOAM/OpenFOAM-6/blob/e05f780ebc875b6fbddff3168e983eb91491828b/src/thermophysicalModels/specie/include/thermoPhysicsTypes.H

#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "constTransport.H"
#include "sutherlandTransport.H"
#include "janafThermo.H"

namespace Foam{



/*
typedef constTransport
    <
        species::thermo
        <
            hConstThermo<perfectGas<specie>>,
            sensibleEnthalpy
        >
    > ThermoType;
*/

typedef
    sutherlandTransport
    <
        species::thermo
        <
            janafThermo
            <
                perfectGas<specie>
            >,
            sensibleEnthalpy
        >
    > ThermoType;



}