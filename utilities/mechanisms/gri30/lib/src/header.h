#ifndef HEAD
#define HEAD
#include <stdlib.h>
#include <math.h>
/** Constant pressure or volume. */
#define CONP
//#define CONV

/** Include mechanism header to get NSP and NN **/
#include "mechanism.h"
// OpenMP
#ifdef _OPENMP
 #include <omp.h>
#else
 #define omp_get_max_threads() 1
 #define omp_get_num_threads() 1
#endif
#endif
