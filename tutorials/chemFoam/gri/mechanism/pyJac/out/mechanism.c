#include "mass_mole.h"
#include <stdio.h>
#include "mechanism.h"
    //apply masking of ICs for cache optimized mechanisms
    void apply_mask(double* y_specs) {
        double temp [NSP];
        memcpy(temp, y_specs, NSP * sizeof(double));
        y_specs[0] = temp[0];
        y_specs[1] = temp[1];
        y_specs[2] = temp[2];
        y_specs[3] = temp[3];
        y_specs[4] = temp[4];
        y_specs[5] = temp[5];
        y_specs[6] = temp[6];
        y_specs[7] = temp[7];
        y_specs[8] = temp[8];
        y_specs[9] = temp[9];
        y_specs[10] = temp[10];
        y_specs[11] = temp[11];
        y_specs[12] = temp[12];
        y_specs[13] = temp[13];
        y_specs[14] = temp[14];
        y_specs[15] = temp[15];
        y_specs[16] = temp[16];
        y_specs[17] = temp[17];
        y_specs[18] = temp[18];
        y_specs[19] = temp[19];
        y_specs[20] = temp[20];
        y_specs[21] = temp[21];
        y_specs[22] = temp[22];
        y_specs[23] = temp[23];
        y_specs[24] = temp[24];
        y_specs[25] = temp[25];
        y_specs[26] = temp[26];
        y_specs[27] = temp[27];
        y_specs[28] = temp[28];
        y_specs[29] = temp[29];
        y_specs[30] = temp[30];
        y_specs[31] = temp[31];
        y_specs[32] = temp[32];
        y_specs[33] = temp[33];
        y_specs[34] = temp[34];
        y_specs[35] = temp[35];
        y_specs[36] = temp[36];
        y_specs[37] = temp[37];
        y_specs[38] = temp[38];
        y_specs[39] = temp[39];
        y_specs[40] = temp[40];
        y_specs[41] = temp[41];
        y_specs[42] = temp[42];
        y_specs[43] = temp[43];
        y_specs[44] = temp[44];
        y_specs[45] = temp[45];
        y_specs[46] = temp[46];
        y_specs[47] = temp[48];
        y_specs[48] = temp[49];
        y_specs[49] = temp[50];
        y_specs[50] = temp[51];
        y_specs[51] = temp[52];
        y_specs[52] = temp[47];
    }
    //reverse masking of ICs for cache optimized mechanisms
    void apply_reverse_mask(double* y_specs) {
        double temp [NSP];
        memcpy(temp, y_specs, NSP * sizeof(double));
        y_specs[0] = temp[0];
        y_specs[1] = temp[1];
        y_specs[2] = temp[2];
        y_specs[3] = temp[3];
        y_specs[4] = temp[4];
        y_specs[5] = temp[5];
        y_specs[6] = temp[6];
        y_specs[7] = temp[7];
        y_specs[8] = temp[8];
        y_specs[9] = temp[9];
        y_specs[10] = temp[10];
        y_specs[11] = temp[11];
        y_specs[12] = temp[12];
        y_specs[13] = temp[13];
        y_specs[14] = temp[14];
        y_specs[15] = temp[15];
        y_specs[16] = temp[16];
        y_specs[17] = temp[17];
        y_specs[18] = temp[18];
        y_specs[19] = temp[19];
        y_specs[20] = temp[20];
        y_specs[21] = temp[21];
        y_specs[22] = temp[22];
        y_specs[23] = temp[23];
        y_specs[24] = temp[24];
        y_specs[25] = temp[25];
        y_specs[26] = temp[26];
        y_specs[27] = temp[27];
        y_specs[28] = temp[28];
        y_specs[29] = temp[29];
        y_specs[30] = temp[30];
        y_specs[31] = temp[31];
        y_specs[32] = temp[32];
        y_specs[33] = temp[33];
        y_specs[34] = temp[34];
        y_specs[35] = temp[35];
        y_specs[36] = temp[36];
        y_specs[37] = temp[37];
        y_specs[38] = temp[38];
        y_specs[39] = temp[39];
        y_specs[40] = temp[40];
        y_specs[41] = temp[41];
        y_specs[42] = temp[42];
        y_specs[43] = temp[43];
        y_specs[44] = temp[44];
        y_specs[45] = temp[45];
        y_specs[46] = temp[46];
        y_specs[47] = temp[52];
        y_specs[48] = temp[47];
        y_specs[49] = temp[48];
        y_specs[50] = temp[49];
        y_specs[51] = temp[50];
        y_specs[52] = temp[51];
    }
void set_same_initial_conditions(int NUM, double** y_host, double** var_host) 
{
    double Xi [NSP] = {0.0};
    //set initial mole fractions here

    //Normalize mole fractions to sum to one
    double Xsum = 0.0;
    for (int j = 0; j < NSP; ++ j) {
        Xsum += Xi[j];
    }
    if (Xsum == 0.0) {
        printf("Use of the set initial conditions function requires user implementation!\n");
        exit(-1);
    }
    for (int j = 0; j < NSP; ++ j) {
        Xi[j] /= Xsum;
    }

    //convert to mass fractions
    double Yi[NSP - 1] = {0.0};
    mole2mass(Xi, Yi);

    //set initial pressure, units [PA]
    double P = 101325.0;
    // set intial temperature, units [K]
    double T0 = 1600;

    (*y_host) = (double*)malloc(NUM * NSP * sizeof(double));
    (*var_host) = (double*)malloc(NUM * sizeof(double));
    //load temperature and mass fractions for all threads (cells)
    for (int i = 0; i < NUM; ++i) {
        (*y_host)[i] = T0;
        //loop through species
        for (int j = 1; j < NSP; ++j) {
            (*y_host)[i + NUM * j] = Yi[j - 1];
        }
    }

#ifdef CONV
    //calculate density
    double rho = getDensity(T0, P, Xi);
#endif

    for (int i = 0; i < NUM; ++i) {
#ifdef CONV
        (*var_host)[i] = rho;
#elif defined(CONP)
        (*var_host)[i] = P;
#endif
    }
}

