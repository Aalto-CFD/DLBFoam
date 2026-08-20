#ifndef MECHANISM_h
#define MECHANISM_h

// The following definitions can be left out due to the prototype function nature.

//Number of species
//#define NSP -1
//Number of variables. NN = NSP + 1 (temperature)
//#define NN -1
//Number of forward reactions
//#define FWD_RATES -1
//Number of reversible reactions
//#define REV_RATES -1
//Number of reactions with pressure modified rates
//#define PRES_MOD_RATES -1

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

//#if defined (RATES_TEST) || defined (PROFILER)
//    void write_jacobian_and_rates_output(int NUM);
//#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

