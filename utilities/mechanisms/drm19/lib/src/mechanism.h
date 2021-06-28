#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 19
/* Species Indexes
0  H2
1  H
2  O
3  O2
4  OH
5  H2O
6  HO2
7  CH2
8  CH2(S)
9  CH3
10  CH4
11  CO
12  CO2
13  HCO
14  CH2O
15  CH3O
16  C2H4
17  C2H5
18  C2H6
19  AR
20  N2
*/

//Number of species
#define NSP 21
//Number of variables. NN = NSP + 1 (temperature)
#define NN 22
//Number of forward reactions
#define FWD_RATES 84
//Number of reversible reactions
#define REV_RATES 84
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 14

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

