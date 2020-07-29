#pragma once
#include <array>
#include <stdlib.h>     //srand, rand
#include "Random.H"

#include "chemistryProblem.H"
#include "chemistrySolution.H"




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

    static constexpr size_t n_species = 53;


    
    chemistryProblem create_heavy_problem() {

        std::vector<double> asd = {
            4.95477e-05, 6.45518e-08, 1.89254e-11, 3.09914e-06, 0.000236918, 6.32301e-05,
            2.47432e-06, 1.68212e-05, 0.00100945,  1.87345e-08, 7.51038e-09, 1.50414e-13,
            3.72389e-12, 1.22597e-09, 3.83324e-11, 6.20063e-07, 6.78764e-05, 0.000456642,
            1.26414e-10, 5.41209e-11, 7.98873e-09, 9.01034e-11, 9.66314e-06, 2.13694e-09,
            1.64621e-06, 1.88428e-09, 2.72789e-07, 3.78312e-09, 2.17603e-07, 7.79603e-09,
            1.21802e-11, 6.30441e-12, 8.01972e-11, 1.6812e-08,  6.16455e-12, 3.69734e-07,
            3.42231e-09, 1.60326e-09, 1.67249e-11, 1.48694e-13, 1.13064e-08, 1.67778e-14,
            9.9342e-15,  1.67998e-09, 6.51119e-12, 5.18643e-09, 4.41158e-12, 0.00475011,
            0,           7.53635e-13, 3.05849e-10, 1.61634e-11, 2.15423e-10};

        chemistryProblem p(n_species);

        std::copy(asd.begin(), asd.end(), p.c.begin());

        p.Ti = 1803.424;
        p.pi = 100000; 
        p.rhoi = 0.1824236;
        p.deltaTChem = 0.000841043;
        p.deltaT = 0.000847458;

        return p;


    }

    chemistryProblem create_light_problem() {

        
        chemistryProblem p(n_species);
        p.c = 0.0;
        p.c[47] = 0.00605864;
        p.Ti = 1985.11;
        p.pi = 99998.7;
        p.rhoi = 0.169723;
        //p.deltaTChem = 1.15399e-05;
        //p.deltaT = 2.98561e-06;
        p.deltaTChem = 0.000841043;
        p.deltaT = 0.000847458;
        return p;


    }
   std::vector<chemistryProblem> create_n_random_heavy(size_t n) {

        std::vector<chemistryProblem> p;

        for (size_t i = 0; i < n;  ++i){
            p.push_back(create_heavy_problem());
        } 
        return p;
    }

     std::vector<chemistryProblem> create_n_random_light(size_t n) {

        std::vector<chemistryProblem> p;

        for (size_t i = 0; i < n;  ++i){
            p.push_back(create_light_problem());
        } 
        return p;
    }

    chemistrySolution create_empty_solution() {
        return chemistrySolution(n_species);
    }

    //static constexpr std::array<size_t, 2> high_species_indices = {4, 47};
    //static constexpr std::array<double, 2> high_species_values = {0.00124406, 0.00475743};

    //static constexpr std::array<double, 2> low_species_values = 

}

void assign_heavy(volScalarField& p, volScalarField& rho, PtrList<volScalarField>& Y,  psiReactionThermo& thermo) {

    auto problem = detail::create_heavy_problem();

    auto field = thermo.W().ref();
    
    forAll(p, celli) {

            p[celli] = problem.pi;
            rho[celli] = problem.rhoi;

            
            for (label i = 0; i < problem.c.size(); i++) { 
                Y[i][celli] = problem.c[i] * field[celli] / problem.rhoi;
            }
            
    }

} 

void assign_light(volScalarField& p, volScalarField& rho, PtrList<volScalarField>& Y,  psiReactionThermo& thermo) {

    auto problem = detail::create_light_problem();

    auto field = thermo.W().ref();
    
    forAll(p, celli) {

            p[celli] = problem.pi;
            rho[celli] = problem.rhoi;

            
            for (label i = 0; i < problem.c.size(); i++) { 
                Y[i][celli] = problem.c[i] * field[celli] / problem.rhoi;
            }
            
    }

} 


struct HighMasterLoadIc{

    
    

    static void set(psiReactionThermo& thermo) {

        /*
        chemistryProblem p;
        if (Pstream::master()){
            p = detail::create_heavy_problem();
            
        }
        else {
            p = detail::create_light_problem();
        }

        //dimensionedScalar(p.dimensions()/dimTime, 0

        //thermo.rho().ref() = dimensionedScalar(thermo.rho().dimension(), p.rhoi);
            //thermo.p().ref() = p.pi;
            //thermo.T().ref() = p.Ti;
        */

    }

};
