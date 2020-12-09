#ifndef CHEM_UTILS_HEAD
#define CHEM_UTILS_HEAD

#include "header.h"

void eval_conc (const double, const double, const double * __restrict__, double * __restrict__, double * __restrict__, double * __restrict__, double * __restrict__);
void eval_conc_rho (const double, const double, const double * __restrict__, double * __restrict__, double * __restrict__, double * __restrict__, double * __restrict__);
void eval_h (const double, double * __restrict__);
void eval_u (const double, double * __restrict__);
void eval_cv (const double, double * __restrict__);
void eval_cp (const double, double * __restrict__);

#endif
