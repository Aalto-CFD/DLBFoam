#ifndef HEAD
#define HEAD

// The header files introduced in this directory work as prototype functions
// for the pyjac C-library, linked in run-time (user-responsibility).

// Due to the function declaration only, no std libs required here.
//#include <stdlib.h>
//#include <math.h>

/** Constant pressure or volume. */
/** For now, pyJacChemistryModel supports only constant pressure mode */
#define CONP
//#define CONV

/** Include mechanism header to get additional NSP and NN **/
#include "mechanism.h"

// OpenMP
#ifdef _OPENMP
 #include <omp.h>
#else
 #define omp_get_max_threads() 1
 #define omp_get_num_threads() 1
#endif
#endif
