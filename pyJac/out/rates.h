#ifndef RATES_HEAD
#define RATES_HEAD

#include "header.h"

void eval_rxn_rates (const double, const double, const double * __restrict__, double * __restrict__, double * __restrict__);
void eval_spec_rates (const double * __restrict__, const double * __restrict__, const double * __restrict__, double * __restrict__, double * __restrict__);
void get_rxn_pres_mod (const double, const double, const double * __restrict__, double * __restrict__);

#endif
