#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 0
/* Species Indexes
0  AR
1  H
2  O
3  OH
4  HO2
5  H2
6  H2O
7  H2O2
8  O2
9  CH2
10  CH2*
11  CH3
12  CH4
13  HCO
14  CH2O
15  CH3O
16  CO
17  CO2
18  C2H2
19  C2H3
20  C2H4
21  C2H5
22  C2H6
23  CH2CHO
24  aC3H5
25  C3H6
26  nC3H7
27  C2H3CHO
28  C4H7
29  C4H81
30  pC4H9
31  C5H9
32  C5H10
33  PXC5H11
34  C6H12
35  PXC6H13
36  C7H14
37  PXC7H15
38  C8H16
39  PXC8H17
40  C9H18
41  PXC9H19
42  C10H20
43  PXC10H21
44  NC12H26
45  PXC12H25
46  SXC12H25
47  S3XC12H25
48  C12H24
49  C12H25O2
50  C12OOH
51  O2C12H24OOH
52  OC12H23OOH
53  N2
*/

//Number of species
#define NSP 54
//Number of variables. NN = NSP + 1 (temperature)
#define NN 55
//Number of forward reactions
#define FWD_RATES 269
//Number of reversible reactions
#define REV_RATES 250
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 31

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

