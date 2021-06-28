#ifndef SPARSE_HEAD
#define SPARSE_HEAD

#define N_A 2809
#include "header.h"

void sparse_multiplier (const double *, const double *, double*);

#ifdef COMPILE_TESTING_METHODS
  int test_sparse_multiplier();
#endif

#endif
