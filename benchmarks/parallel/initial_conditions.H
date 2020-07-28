#pragma once


#include "Random.H"

struct random_rho_T_Y_initial_condition{

    template<class State>
    static void set(State& g_state) {


        Random  randomizer1(541); 
        Random  randomizer2(44); 
        Random  randomizer3(10);

        forAll(g_state.rho, i) {
            g_state.rho[i] = 1.2 + randomizer1.sample01<scalar>() * 3;
            g_state.thermo.T()[i] = 500 + randomizer2.sample01<scalar>() * 100;
        }

        forAll(g_state.rho, celli) {

            for (label i = 0; i < g_state.Y.size(); i++) { 
                g_state.Y[i][celli] = randomizer3.sample01<scalar>();
                //c_[i] = Y_[i][celli]; 
            }

        }

    }

};
