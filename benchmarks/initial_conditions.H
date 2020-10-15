#pragma once
#include <array>
#include <stdlib.h>     //srand, rand
#include "Random.H"

#include "ChemistryProblem.H"
#include "ChemistrySolution.H"




/*
scalarField c;
    scalar      Ti;
    scalar      pi;
    scalar      rhoi;
    scalar      deltaTChem;
    scalar      deltaT;
    scalar      cpuTime;
    label       cellid;

*/


namespace detail{



    static double random_double(double min, double max) {
        double f = double(rand()) / RAND_MAX;
        return min + f * (max - min);
    }

    static constexpr int n_species = 53;
    static constexpr double FLOWTIMESTEP = 1E-7;
    static constexpr double HEAVY_Y_VALUE = 0.001;
    static constexpr double LIGHT_Y_VALUE = 0.0001;


    ChemistryProblem create_heavy_problem() {


        ChemistryProblem p(n_species);
        p.c = 0.1;

        p.Ti = 1800;
        p.pi = 100000;
        p.rhoi = 0.1824236;
        p.deltaTChem = 0.000841043;
        p.deltaT = 0.000847458;

        return p;


    }

    ChemistryProblem create_light_problem() {

        auto p = create_heavy_problem();
        p.c = 0.0001;
        return p;


    }
   std::vector<ChemistryProblem> create_n_random_heavy(size_t n) {

        std::vector<ChemistryProblem> p;

        for (size_t i = 0; i < n;  ++i){
            p.push_back(create_heavy_problem());
        }
        return p;
    }

     std::vector<ChemistryProblem> create_n_random_light(size_t n) {

        std::vector<ChemistryProblem> p;

        for (size_t i = 0; i < n;  ++i){
            p.push_back(create_light_problem());
        }
        return p;
    }

    ChemistrySolution create_empty_solution() {
        return ChemistrySolution(n_species);
    }

void generic_assign_all(volScalarField& rho, psiReactionThermo& thermo, scalar Tv, scalar pv, scalar rhov, scalar Yv) {


    PtrList<volScalarField>& Y = thermo.composition().Y();

    volScalarField& T = thermo.T();
    volScalarField& p = thermo.p();
    forAll(p, celli) {

            p[celli] = pv;
            rho[celli] = rhov;
            T[celli] = Tv;

            for (label i = 0; i < n_species; i++) {
                Y[i][celli] = Yv; //detail::random_double(0.05, 0.1);
            }

    }


}

template<class fractionType>
void generic_assign_all_random(volScalarField& rho, psiReactionThermo& thermo, scalar Tv, scalar pv, scalar rhov, const fractionType& Yv) {


    PtrList<volScalarField>& Y = thermo.composition().Y();

    volScalarField& T = thermo.T();
    volScalarField& p = thermo.p();
    forAll(p, celli) {

            p[celli] = pv;
            rho[celli] = rhov;
            T[celli] = Tv;

            for (label i = 0; i < n_species; i++) {
                Y[i][celli] = Yv(); //detail::random_double(0.05, 0.1);
            }

    }


}

}




void set_all_heavy(volScalarField& rho,  psiReactionThermo& thermo) {

    auto problem = detail::create_heavy_problem();

    detail::generic_assign_all(rho, thermo, problem.Ti, problem.pi, problem.rhoi, detail::HEAVY_Y_VALUE);

}


void set_all_light( volScalarField& rho,  psiReactionThermo& thermo) {

    auto problem = detail::create_light_problem();

    detail::generic_assign_all(rho, thermo, problem.Ti, problem.pi, problem.rhoi, detail::LIGHT_Y_VALUE);

}



void set_master_heavy( volScalarField& rho,  psiReactionThermo& thermo) {


    if (Pstream::master()){
        set_all_heavy(rho, thermo);
    }

    else {
        set_all_light(rho, thermo);
    }
}

void set_every_n_heavy( volScalarField& rho,  psiReactionThermo& thermo, size_t n) {




    if (Pstream::myProcNo() % n == 0){
        set_all_heavy(rho, thermo);
    }

    else {
        set_all_light(rho, thermo);
    }

}


void set_every_n_random_heavy( volScalarField& rho,  psiReactionThermo& thermo, size_t n) {

    auto problem = detail::create_light_problem();

    struct randomizer_heavy {

        scalar operator()() const{
            return detail::random_double(detail::HEAVY_Y_VALUE*0.5, detail::HEAVY_Y_VALUE*1.5);
        }
    };

    struct randomizer_light {

        scalar operator()() const{
            return detail::random_double(detail::LIGHT_Y_VALUE*0.5, detail::LIGHT_Y_VALUE*1.5);
        }
    };

    if (Pstream::myProcNo() % n == 0){
        detail::generic_assign_all_random(rho, thermo, problem.Ti, problem.pi, problem.rhoi, randomizer_heavy());
    }

    else {

        detail::generic_assign_all_random(rho, thermo, problem.Ti, problem.pi, problem.rhoi, randomizer_light());

    }

}


void set_all_random( volScalarField& rho,  psiReactionThermo& thermo) {

    auto problem = detail::create_light_problem();

    struct randomizer {

        scalar operator()() const{
            return detail::random_double(detail::LIGHT_Y_VALUE, detail::HEAVY_Y_VALUE);
        }
    };



    detail::generic_assign_all_random(rho, thermo, problem.Ti, problem.pi, problem.rhoi, randomizer());




}