#include "header.h"
#include "chem_utils.h"
#include "rates.h"

#if defined(CONP)

void dydt (const double t, const double pres, const double * __restrict__ y, double * __restrict__ dy) {

  // species molar concentrations
  double conc[21];
  double y_N;
  double mw_avg;
  double rho;
  eval_conc (y[0], pres, &y[1], &y_N, &mw_avg, &rho, conc);

  // local arrays holding reaction rates
  double fwd_rates[84];
  double rev_rates[84];
  eval_rxn_rates (y[0], pres, conc, fwd_rates, rev_rates);

  // get pressure modifications to reaction rates
  double pres_mod[14];
  get_rxn_pres_mod (y[0], pres, conc, pres_mod);

  // evaluate species molar net production rates
  double dy_N;
  eval_spec_rates (fwd_rates, rev_rates, pres_mod, &dy[1], &dy_N);
  // local array holding constant pressure specific heat
  double cp[21];
  eval_cp (y[0], cp);

  // constant pressure mass-average specific heat
  double cp_avg = (cp[0] * y[1]) + (cp[1] * y[2]) + (cp[2] * y[3]) + (cp[3] * y[4])
              + (cp[4] * y[5]) + (cp[5] * y[6]) + (cp[6] * y[7]) + (cp[7] * y[8])
              + (cp[8] * y[9]) + (cp[9] * y[10]) + (cp[10] * y[11]) + (cp[11] * y[12])
              + (cp[12] * y[13]) + (cp[13] * y[14]) + (cp[14] * y[15]) + (cp[15] * y[16])
              + (cp[16] * y[17]) + (cp[17] * y[18]) + (cp[18] * y[19]) + (cp[19] * y[20]) + (cp[20] * y_N);

  // local array for species enthalpies
  double h[21];
  eval_h(y[0], h);
  // rate of change of temperature
  dy[0] = (-1.0 / (rho * cp_avg)) * ((dy[1] * h[0] * 2.0158800000000001e+00)
        + (dy[2] * h[1] * 1.0079400000000001e+00) + (dy[3] * h[2] * 1.5999400000000000e+01)
        + (dy[4] * h[3] * 3.1998799999999999e+01) + (dy[5] * h[4] * 1.7007339999999999e+01)
        + (dy[6] * h[5] * 1.8015280000000001e+01) + (dy[7] * h[6] * 3.3006740000000001e+01)
        + (dy[8] * h[7] * 1.4026879999999998e+01) + (dy[9] * h[8] * 1.4026879999999998e+01)
        + (dy[10] * h[9] * 1.5034820000000000e+01) + (dy[11] * h[10] * 1.6042760000000001e+01)
        + (dy[12] * h[11] * 2.8010399999999997e+01) + (dy[13] * h[12] * 4.4009799999999998e+01)
        + (dy[14] * h[13] * 2.9018339999999998e+01) + (dy[15] * h[14] * 3.0026280000000000e+01)
        + (dy[16] * h[15] * 3.1034219999999998e+01) + (dy[17] * h[16] * 2.8053759999999997e+01)
        + (dy[18] * h[17] * 2.9061699999999998e+01) + (dy[19] * h[18] * 3.0069640000000000e+01));

  // calculate rate of change of species mass fractions
  dy[1] *= (2.0158800000000001e+00 / rho);
  dy[2] *= (1.0079400000000001e+00 / rho);
  dy[3] *= (1.5999400000000000e+01 / rho);
  dy[4] *= (3.1998799999999999e+01 / rho);
  dy[5] *= (1.7007339999999999e+01 / rho);
  dy[6] *= (1.8015280000000001e+01 / rho);
  dy[7] *= (3.3006740000000001e+01 / rho);
  dy[8] *= (1.4026879999999998e+01 / rho);
  dy[9] *= (1.4026879999999998e+01 / rho);
  dy[10] *= (1.5034820000000000e+01 / rho);
  dy[11] *= (1.6042760000000001e+01 / rho);
  dy[12] *= (2.8010399999999997e+01 / rho);
  dy[13] *= (4.4009799999999998e+01 / rho);
  dy[14] *= (2.9018339999999998e+01 / rho);
  dy[15] *= (3.0026280000000000e+01 / rho);
  dy[16] *= (3.1034219999999998e+01 / rho);
  dy[17] *= (2.8053759999999997e+01 / rho);
  dy[18] *= (2.9061699999999998e+01 / rho);
  dy[19] *= (3.0069640000000000e+01 / rho);
  dy[20] *= (3.9948000000000000e+01 / rho);

} // end dydt

#elif defined(CONV)

void dydt (const double t, const double rho, const double * __restrict__ y, double * __restrict__ dy) {

  // species molar concentrations
  double conc[21];
  double y_N;
  double mw_avg;
  double pres;
  eval_conc_rho (y[0]rho, &y[1], &y_N, &mw_avg, &pres, conc);

  // local arrays holding reaction rates
  double fwd_rates[84];
  double rev_rates[84];
  eval_rxn_rates (y[0], pres, conc, fwd_rates, rev_rates);

  // get pressure modifications to reaction rates
  double pres_mod[14];
  get_rxn_pres_mod (y[0], pres, conc, pres_mod);

  // evaluate species molar net production rates
  double dy_N;  eval_spec_rates (fwd_rates, rev_rates, pres_mod, &dy[1], &dy_N);

  double cv[21];
  eval_cv(y[0], cv);

  // constant volume mass-average specific heat
  double cv_avg = (cv[0] * y[1]) + (cv[1] * y[2]) + (cv[2] * y[3]) + (cv[3] * y[4])
              + (cv[4] * y[5]) + (cv[5] * y[6]) + (cv[6] * y[7]) + (cv[7] * y[8])
              + (cv[8] * y[9]) + (cv[9] * y[10]) + (cv[10] * y[11]) + (cv[11] * y[12])
              + (cv[12] * y[13]) + (cv[13] * y[14]) + (cv[14] * y[15]) + (cv[15] * y[16])
              + (cv[16] * y[17]) + (cv[17] * y[18]) + (cv[18] * y[19]) + (cv[19] * y[20])(cv[20] * y_N);

  // local array for species internal energies
  double u[21];
  eval_u (y[0], u);

  // rate of change of temperature
  dy[0] = (-1.0 / (rho * cv_avg)) * ((dy[1] * u[0] * 2.0158800000000001e+00)
        + (dy[2] * u[1] * 1.0079400000000001e+00) + (dy[3] * u[2] * 1.5999400000000000e+01)
        + (dy[4] * u[3] * 3.1998799999999999e+01) + (dy[5] * u[4] * 1.7007339999999999e+01)
        + (dy[6] * u[5] * 1.8015280000000001e+01) + (dy[7] * u[6] * 3.3006740000000001e+01)
        + (dy[8] * u[7] * 1.4026879999999998e+01) + (dy[9] * u[8] * 1.4026879999999998e+01)
        + (dy[10] * u[9] * 1.5034820000000000e+01) + (dy[11] * u[10] * 1.6042760000000001e+01)
        + (dy[12] * u[11] * 2.8010399999999997e+01) + (dy[13] * u[12] * 4.4009799999999998e+01)
        + (dy[14] * u[13] * 2.9018339999999998e+01) + (dy[15] * u[14] * 3.0026280000000000e+01)
        + (dy[16] * u[15] * 3.1034219999999998e+01) + (dy[17] * u[16] * 2.8053759999999997e+01)
        + (dy[18] * u[17] * 2.9061699999999998e+01) + (dy[19] * u[18] * 3.0069640000000000e+01));

  // calculate rate of change of species mass fractions
  dy[1] *= (2.0158800000000001e+00 / rho);
  dy[2] *= (1.0079400000000001e+00 / rho);
  dy[3] *= (1.5999400000000000e+01 / rho);
  dy[4] *= (3.1998799999999999e+01 / rho);
  dy[5] *= (1.7007339999999999e+01 / rho);
  dy[6] *= (1.8015280000000001e+01 / rho);
  dy[7] *= (3.3006740000000001e+01 / rho);
  dy[8] *= (1.4026879999999998e+01 / rho);
  dy[9] *= (1.4026879999999998e+01 / rho);
  dy[10] *= (1.5034820000000000e+01 / rho);
  dy[11] *= (1.6042760000000001e+01 / rho);
  dy[12] *= (2.8010399999999997e+01 / rho);
  dy[13] *= (4.4009799999999998e+01 / rho);
  dy[14] *= (2.9018339999999998e+01 / rho);
  dy[15] *= (3.0026280000000000e+01 / rho);
  dy[16] *= (3.1034219999999998e+01 / rho);
  dy[17] *= (2.8053759999999997e+01 / rho);
  dy[18] *= (2.9061699999999998e+01 / rho);
  dy[19] *= (3.0069640000000000e+01 / rho);
  dy[20] *= (3.9948000000000000e+01 / rho);

} // end dydt

#endif
