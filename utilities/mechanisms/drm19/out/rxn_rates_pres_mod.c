#include <math.h>
#include "header.h"
#include "rates.h"

void get_rxn_pres_mod (const double T, const double pres, const double * __restrict__ C, double * __restrict__ pres_mod) {
  // third body variable declaration
  double thd;

  // pressure dependence variable declarations
  double k0;
  double kinf;
  double Pr;

  // troe variable declarations
  double logFcent;
  double A;
  double B;

  double logT = log(T);
  double m = pres / (8.31446210e+03 * T);

  // reaction 0;
  pres_mod[0] = m + 1.0 * C[0] + 5.0 * C[5] + 1.0 * C[10] + 0.5 * C[11] + 1.0 * C[12] + 2.0 * C[18] - 0.3 * C[19];

  // reaction 7;
  pres_mod[1] = m + 1.0 * C[0] + 5.0 * C[3] + 5.0 * C[5] + 1.0 * C[10] + 0.5 * C[11] + 2.5 * C[12] + 2.0 * C[18] - 0.5 * C[19];

  // reaction 16;
  pres_mod[2] = m - 1.0 * C[3] - 1.0 * C[5] - 0.25 * C[11] + 0.5 * C[12] + 0.5 * C[18] - 1.0 * C[20] - 1.0 * C[19];

  // reaction 22;
  pres_mod[3] = m - 1.0 * C[0] - 1.0 * C[5] + 1.0 * C[10] - 1.0 * C[12] + 2.0 * C[18] - 0.37 * C[19];

  // reaction 26;
  pres_mod[4] = m - 0.27 * C[0] + 2.65 * C[5] + 1.0 * C[10] + 2.0 * C[18] - 0.62 * C[19];

  // reaction 29;
  thd = m + 1.0 * C[0] + 5.0 * C[5] + 1.0 * C[10] + 0.5 * C[11] + 1.0 * C[12] + 2.0 * C[18] - 0.3 * C[19];
  k0 = exp(4.9517437762680643e+01 - 3.14 * logT - (6.1896006477677008e+02 / T));
  kinf = exp(3.0849896940796750e+01 - 0.8 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(3.20000000e-01 * exp(-T / 7.80000000e+01) + 6.80000000e-01 * exp(-T / 1.99500000e+03) + exp(-5.59000000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[5] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 30;
  thd = m + 1.0 * C[0] + 5.0 * C[5] + 1.0 * C[10] + 0.5 * C[11] + 1.0 * C[12] + 2.0 * C[18] - 0.3 * C[19];
  k0 = exp(6.3076845661346454e+01 - 4.76 * logT - (1.2278557382563567e+03 / T));
  kinf = exp(3.0172623109393093e+01 - 0.63 * logT - (1.9273309334105929e+02 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.17000000e-01 * exp(-T / 7.40000000e+01) + 7.83000000e-01 * exp(-T / 2.94100000e+03) + exp(-6.96400000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[6] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 32;
  thd = m + 1.0 * C[0] + 5.0 * C[5] + 1.0 * C[10] + 0.5 * C[11] + 1.0 * C[12] + 2.0 * C[18] - 0.3 * C[19];
  k0 = exp(4.1746636266343160e+01 - 2.57 * logT - (7.1708787992430678e+02 / T));
  kinf = exp(2.0809443533187462e+01 + 0.48 * logT - (-1.3083708686338230e+02 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.17600000e-01 * exp(-T / 2.71000000e+02) + 7.82400000e-01 * exp(-T / 2.75500000e+03) + exp(-6.57000000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[7] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 34;
  thd = m + 1.0 * C[0] + 5.0 * C[5] + 1.0 * C[10] + 0.5 * C[11] + 1.0 * C[12] + 2.0 * C[18];
  k0 = exp(5.6050499592221364e+01 - 4.8 * logT - (2.7979007806169443e+03 / T));
  kinf = exp(2.0107079697522593e+01 + 0.454 * logT - (1.3083708686338230e+03 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.42000000e-01 * exp(-T / 9.40000000e+01) + 7.58000000e-01 * exp(-T / 1.55500000e+03) + exp(-4.20000000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[8] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 37;
  thd = m + 1.0 * C[0] + 5.0 * C[5] + 1.0 * C[10] + 0.5 * C[11] + 1.0 * C[12] + 2.0 * C[18] - 0.3 * C[19];
  k0 = exp(8.3075384904579593e+01 - 7.62 * logT - (3.5074403670683637e+03 / T));
  kinf = exp(2.0800226878082540e+01 + 0.454 * logT - (9.1585960804367608e+02 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.47000000e-02 * exp(-T / 2.10000000e+02) + 9.75300000e-01 * exp(-T / 9.84000000e+02) + exp(-4.37400000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[9] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 38;
  thd = m + 1.0 * C[0] + 5.0 * C[5] + 1.0 * C[10] + 0.5 * C[11] + 1.0 * C[12] + 2.0 * C[18] - 0.3 * C[19];
  k0 = exp(8.1278612893528006e+01 - 7.08 * logT - (3.3640227910835024e+03 / T));
  kinf = exp(3.3886771157681913e+01 - 0.99 * logT - (7.9508691247747697e+02 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(1.57800000e-01 * exp(-T / 1.25000000e+02) + 8.42200000e-01 * exp(-T / 2.21900000e+03) + exp(-6.88200000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[10] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 40;
  thd = m + 1.0 * C[0] + 5.0 * C[5] + 1.0 * C[10] + 0.5 * C[11] + 1.0 * C[12] + 2.0 * C[18] - 0.3 * C[19];
  k0 = exp(4.9977627770478051e+01 - 3.42 * logT - (4.2446570295870370e+04 / T));
  kinf = exp(1.0668955394675699e+01 + 1.5 * logT - (4.0056277362789348e+04 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(6.80000000e-02 * exp(-T / 1.97000000e+02) + 9.32000000e-01 * exp(-T / 1.54000000e+03) + exp(-1.03000000e+04 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[11] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 74;
  thd = m + 1.0 * C[0] + 5.0 * C[5] + 1.0 * C[10] + 0.5 * C[11] + 1.0 * C[12] + 2.0 * C[18] - 0.3 * C[19];
  k0 = exp(1.0188472363832375e+02 - 9.67 * logT - (3.1300256934239915e+03 / T));
  kinf = exp(3.0685022297606515e+01 - 0.97 * logT - (3.1199613021268084e+02 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(4.67500000e-01 * exp(-T / 1.51000000e+02) + 5.32500000e-01 * exp(-T / 1.03800000e+03) + exp(-4.97000000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[12] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 80;
  pres_mod[13] = m + 1.0 * C[0] - 1.0 * C[5] + 1.0 * C[10] + 0.5 * C[11] + 1.0 * C[12] + 2.0 * C[18];

} // end get_rxn_pres_mod

