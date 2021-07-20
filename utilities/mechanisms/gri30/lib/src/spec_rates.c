#include "header.h"
#include "rates.h"

void eval_spec_rates (const double * __restrict__ fwd_rates, const double * __restrict__ rev_rates, const double * __restrict__ pres_mod, double * __restrict__ sp_rates, double * __restrict__ dy_N) {
  //rxn 0
  //sp 2
  sp_rates[2] = -2.0 * (fwd_rates[0] - rev_rates[0]) * pres_mod[0];
  //sp 3
  sp_rates[3] = (fwd_rates[0] - rev_rates[0]) * pres_mod[0];

  //rxn 1
  //sp 1
  sp_rates[1] = -(fwd_rates[1] - rev_rates[1]) * pres_mod[1];
  //sp 2
  sp_rates[2] -= (fwd_rates[1] - rev_rates[1]) * pres_mod[1];
  //sp 4
  sp_rates[4] = (fwd_rates[1] - rev_rates[1]) * pres_mod[1];

  //rxn 2
  //sp 0
  sp_rates[0] = -(fwd_rates[2] - rev_rates[2]);
  //sp 1
  sp_rates[1] += (fwd_rates[2] - rev_rates[2]);
  //sp 2
  sp_rates[2] -= (fwd_rates[2] - rev_rates[2]);
  //sp 4
  sp_rates[4] += (fwd_rates[2] - rev_rates[2]);

  //rxn 3
  //sp 2
  sp_rates[2] -= (fwd_rates[3] - rev_rates[3]);
  //sp 3
  sp_rates[3] += (fwd_rates[3] - rev_rates[3]);
  //sp 4
  sp_rates[4] += (fwd_rates[3] - rev_rates[3]);
  //sp 6
  sp_rates[6] = -(fwd_rates[3] - rev_rates[3]);

  //rxn 4
  //sp 2
  sp_rates[2] -= (fwd_rates[4] - rev_rates[4]);
  //sp 4
  sp_rates[4] += (fwd_rates[4] - rev_rates[4]);
  //sp 6
  sp_rates[6] += (fwd_rates[4] - rev_rates[4]);
  //sp 7
  sp_rates[7] = -(fwd_rates[4] - rev_rates[4]);

  //rxn 5
  //sp 9
  sp_rates[9] = -(fwd_rates[5] - rev_rates[5]);
  //sp 2
  sp_rates[2] -= (fwd_rates[5] - rev_rates[5]);
  //sp 14
  sp_rates[14] = (fwd_rates[5] - rev_rates[5]);
  //sp 1
  sp_rates[1] += (fwd_rates[5] - rev_rates[5]);

  //rxn 6
  //sp 16
  sp_rates[16] = (fwd_rates[6] - rev_rates[6]);
  //sp 1
  sp_rates[1] += (fwd_rates[6] - rev_rates[6]);
  //sp 10
  sp_rates[10] = -(fwd_rates[6] - rev_rates[6]);
  //sp 2
  sp_rates[2] -= (fwd_rates[6] - rev_rates[6]);

  //rxn 7
  //sp 0
  sp_rates[0] += (fwd_rates[7] - rev_rates[7]);
  //sp 2
  sp_rates[2] -= (fwd_rates[7] - rev_rates[7]);
  //sp 11
  sp_rates[11] = -(fwd_rates[7] - rev_rates[7]);
  //sp 14
  sp_rates[14] += (fwd_rates[7] - rev_rates[7]);

  //rxn 8
  //sp 16
  sp_rates[16] += (fwd_rates[8] - rev_rates[8]);
  //sp 1
  sp_rates[1] += (fwd_rates[8] - rev_rates[8]);
  //sp 2
  sp_rates[2] -= (fwd_rates[8] - rev_rates[8]);
  //sp 11
  sp_rates[11] -= (fwd_rates[8] - rev_rates[8]);

  //rxn 9
  //sp 17
  sp_rates[17] = (fwd_rates[9] - rev_rates[9]);
  //sp 2
  sp_rates[2] -= (fwd_rates[9] - rev_rates[9]);
  //sp 12
  sp_rates[12] = -(fwd_rates[9] - rev_rates[9]);
  //sp 1
  sp_rates[1] += (fwd_rates[9] - rev_rates[9]);

  //rxn 10
  //sp 4
  sp_rates[4] += (fwd_rates[10] - rev_rates[10]);
  //sp 2
  sp_rates[2] -= (fwd_rates[10] - rev_rates[10]);
  //sp 12
  sp_rates[12] += (fwd_rates[10] - rev_rates[10]);
  //sp 13
  sp_rates[13] = -(fwd_rates[10] - rev_rates[10]);

  //rxn 11
  //sp 2
  sp_rates[2] -= (fwd_rates[11] - rev_rates[11]) * pres_mod[2];
  //sp 14
  sp_rates[14] -= (fwd_rates[11] - rev_rates[11]) * pres_mod[2];
  //sp 15
  sp_rates[15] = (fwd_rates[11] - rev_rates[11]) * pres_mod[2];

  //rxn 12
  //sp 16
  sp_rates[16] -= (fwd_rates[12] - rev_rates[12]);
  //sp 2
  sp_rates[2] -= (fwd_rates[12] - rev_rates[12]);
  //sp 4
  sp_rates[4] += (fwd_rates[12] - rev_rates[12]);
  //sp 14
  sp_rates[14] += (fwd_rates[12] - rev_rates[12]);

  //rxn 13
  //sp 16
  sp_rates[16] -= (fwd_rates[13] - rev_rates[13]);
  //sp 1
  sp_rates[1] += (fwd_rates[13] - rev_rates[13]);
  //sp 2
  sp_rates[2] -= (fwd_rates[13] - rev_rates[13]);
  //sp 15
  sp_rates[15] += (fwd_rates[13] - rev_rates[13]);

  //rxn 14
  //sp 16
  sp_rates[16] += (fwd_rates[14] - rev_rates[14]);
  //sp 17
  sp_rates[17] -= (fwd_rates[14] - rev_rates[14]);
  //sp 2
  sp_rates[2] -= (fwd_rates[14] - rev_rates[14]);
  //sp 4
  sp_rates[4] += (fwd_rates[14] - rev_rates[14]);

  //rxn 15
  //sp 17
  sp_rates[17] += (fwd_rates[15] - rev_rates[15]);
  //sp 18
  sp_rates[18] = -(fwd_rates[15] - rev_rates[15]);
  //sp 4
  sp_rates[4] += (fwd_rates[15] - rev_rates[15]);
  //sp 2
  sp_rates[2] -= (fwd_rates[15] - rev_rates[15]);

  //rxn 16
  //sp 17
  sp_rates[17] += (fwd_rates[16] - rev_rates[16]);
  //sp 2
  sp_rates[2] -= (fwd_rates[16] - rev_rates[16]);
  //sp 19
  sp_rates[19] = -(fwd_rates[16] - rev_rates[16]);
  //sp 4
  sp_rates[4] += (fwd_rates[16] - rev_rates[16]);

  //rxn 17
  //sp 4
  sp_rates[4] += (fwd_rates[17] - rev_rates[17]);
  //sp 2
  sp_rates[2] -= (fwd_rates[17] - rev_rates[17]);
  //sp 20
  sp_rates[20] = -(fwd_rates[17] - rev_rates[17]);
  //sp 18
  sp_rates[18] += (fwd_rates[17] - rev_rates[17]);

  //rxn 18
  //sp 4
  sp_rates[4] += (fwd_rates[18] - rev_rates[18]);
  //sp 2
  sp_rates[2] -= (fwd_rates[18] - rev_rates[18]);
  //sp 19
  sp_rates[19] += (fwd_rates[18] - rev_rates[18]);
  //sp 20
  sp_rates[20] -= (fwd_rates[18] - rev_rates[18]);

  //rxn 19
  //sp 9
  sp_rates[9] += (fwd_rates[19] - rev_rates[19]);
  //sp 2
  sp_rates[2] -= (fwd_rates[19] - rev_rates[19]);
  //sp 21
  sp_rates[21] = -(fwd_rates[19] - rev_rates[19]);
  //sp 14
  sp_rates[14] += (fwd_rates[19] - rev_rates[19]);

  //rxn 20
  //sp 1
  sp_rates[1] += (fwd_rates[20] - rev_rates[20]);
  //sp 2
  sp_rates[2] -= (fwd_rates[20] - rev_rates[20]);
  //sp 27
  sp_rates[27] = (fwd_rates[20] - rev_rates[20]);
  //sp 22
  sp_rates[22] = -(fwd_rates[20] - rev_rates[20]);

  //rxn 21
  //sp 2
  sp_rates[2] -= (fwd_rates[21] - rev_rates[21]);
  //sp 4
  sp_rates[4] += (fwd_rates[21] - rev_rates[21]);
  //sp 21
  sp_rates[21] += (fwd_rates[21] - rev_rates[21]);
  //sp 22
  sp_rates[22] -= (fwd_rates[21] - rev_rates[21]);

  //rxn 22
  //sp 2
  sp_rates[2] -= (fwd_rates[22] - rev_rates[22]);
  //sp 14
  sp_rates[14] += (fwd_rates[22] - rev_rates[22]);
  //sp 10
  sp_rates[10] += (fwd_rates[22] - rev_rates[22]);
  //sp 22
  sp_rates[22] -= (fwd_rates[22] - rev_rates[22]);

  //rxn 23
  //sp 1
  sp_rates[1] += (fwd_rates[23] - rev_rates[23]);
  //sp 2
  sp_rates[2] -= (fwd_rates[23] - rev_rates[23]);
  //sp 28
  sp_rates[28] = (fwd_rates[23] - rev_rates[23]);
  //sp 23
  sp_rates[23] = -(fwd_rates[23] - rev_rates[23]);

  //rxn 24
  //sp 24
  sp_rates[24] = -(fwd_rates[24] - rev_rates[24]);
  //sp 16
  sp_rates[16] += (fwd_rates[24] - rev_rates[24]);
  //sp 2
  sp_rates[2] -= (fwd_rates[24] - rev_rates[24]);
  //sp 12
  sp_rates[12] += (fwd_rates[24] - rev_rates[24]);

  //rxn 25
  //sp 25
  sp_rates[25] = -(fwd_rates[25] - rev_rates[25]);
  //sp 2
  sp_rates[2] -= (fwd_rates[25] - rev_rates[25]);
  //sp 12
  sp_rates[12] += (fwd_rates[25] - rev_rates[25]);
  //sp 17
  sp_rates[17] += (fwd_rates[25] - rev_rates[25]);

  //rxn 26
  //sp 25
  sp_rates[25] += (fwd_rates[26] - rev_rates[26]);
  //sp 26
  sp_rates[26] = -(fwd_rates[26] - rev_rates[26]);
  //sp 4
  sp_rates[4] += (fwd_rates[26] - rev_rates[26]);
  //sp 2
  sp_rates[2] -= (fwd_rates[26] - rev_rates[26]);

  //rxn 27
  //sp 1
  sp_rates[1] += (fwd_rates[27] - rev_rates[27]);
  //sp 2
  sp_rates[2] -= (fwd_rates[27] - rev_rates[27]);
  //sp 27
  sp_rates[27] -= (fwd_rates[27] - rev_rates[27]);
  //sp 14
  sp_rates[14] += 2.0 * (fwd_rates[27] - rev_rates[27]);

  //rxn 28
  //sp 4
  sp_rates[4] += (fwd_rates[28] - rev_rates[28]);
  //sp 2
  sp_rates[2] -= (fwd_rates[28] - rev_rates[28]);
  //sp 27
  sp_rates[27] += (fwd_rates[28] - rev_rates[28]);
  //sp 28
  sp_rates[28] -= (fwd_rates[28] - rev_rates[28]);

  //rxn 29
  //sp 2
  sp_rates[2] -= (fwd_rates[29] - rev_rates[29]);
  //sp 28
  sp_rates[28] -= (fwd_rates[29] - rev_rates[29]);
  //sp 10
  sp_rates[10] += (fwd_rates[29] - rev_rates[29]);
  //sp 15
  sp_rates[15] += (fwd_rates[29] - rev_rates[29]);

  //rxn 30
  //sp 2
  sp_rates[2] += (fwd_rates[30] - rev_rates[30]);
  //sp 3
  sp_rates[3] -= (fwd_rates[30] - rev_rates[30]);
  //sp 14
  sp_rates[14] -= (fwd_rates[30] - rev_rates[30]);
  //sp 15
  sp_rates[15] += (fwd_rates[30] - rev_rates[30]);

  //rxn 31
  //sp 16
  sp_rates[16] += (fwd_rates[31] - rev_rates[31]);
  //sp 17
  sp_rates[17] -= (fwd_rates[31] - rev_rates[31]);
  //sp 3
  sp_rates[3] -= (fwd_rates[31] - rev_rates[31]);
  //sp 6
  sp_rates[6] += (fwd_rates[31] - rev_rates[31]);

  //rxn 32
  //sp 1
  sp_rates[1] -= (fwd_rates[32] - rev_rates[32]) * pres_mod[3];
  //sp 3
  sp_rates[3] -= (fwd_rates[32] - rev_rates[32]) * pres_mod[3];
  //sp 6
  sp_rates[6] += (fwd_rates[32] - rev_rates[32]) * pres_mod[3];

  //rxn 33
  //sp 1
  sp_rates[1] -= (fwd_rates[33] - rev_rates[33]);
  //sp 3
  sp_rates[3] -= (fwd_rates[33] - rev_rates[33]);
  //sp 6
  sp_rates[6] += (fwd_rates[33] - rev_rates[33]);

  //rxn 34
  //sp 1
  sp_rates[1] -= (fwd_rates[34] - rev_rates[34]);
  //sp 3
  sp_rates[3] -= (fwd_rates[34] - rev_rates[34]);
  //sp 6
  sp_rates[6] += (fwd_rates[34] - rev_rates[34]);

  //rxn 35
  //sp 1
  sp_rates[1] -= (fwd_rates[35] - rev_rates[35]);
  //sp 3
  sp_rates[3] -= (fwd_rates[35] - rev_rates[35]);
  //sp 6
  sp_rates[6] += (fwd_rates[35] - rev_rates[35]);

  //rxn 36
  //sp 1
  sp_rates[1] -= (fwd_rates[36] - rev_rates[36]);
  //sp 3
  sp_rates[3] -= (fwd_rates[36] - rev_rates[36]);
  //sp 6
  sp_rates[6] += (fwd_rates[36] - rev_rates[36]);

  //rxn 37
  //sp 1
  sp_rates[1] -= (fwd_rates[37] - rev_rates[37]);
  //sp 2
  sp_rates[2] += (fwd_rates[37] - rev_rates[37]);
  //sp 3
  sp_rates[3] -= (fwd_rates[37] - rev_rates[37]);
  //sp 4
  sp_rates[4] += (fwd_rates[37] - rev_rates[37]);

  //rxn 38
  //sp 0
  sp_rates[0] += (fwd_rates[38] - rev_rates[38]) * pres_mod[4];
  //sp 1
  sp_rates[1] -= 2.0 * (fwd_rates[38] - rev_rates[38]) * pres_mod[4];

  //rxn 39
  //sp 0
  sp_rates[0] += (fwd_rates[39] - rev_rates[39]);
  //sp 1
  sp_rates[1] -= 2.0 * (fwd_rates[39] - rev_rates[39]);

  //rxn 40
  //sp 0
  sp_rates[0] += (fwd_rates[40] - rev_rates[40]);
  //sp 1
  sp_rates[1] -= 2.0 * (fwd_rates[40] - rev_rates[40]);

  //rxn 41
  //sp 0
  sp_rates[0] += (fwd_rates[41] - rev_rates[41]);
  //sp 1
  sp_rates[1] -= 2.0 * (fwd_rates[41] - rev_rates[41]);

  //rxn 42
  //sp 1
  sp_rates[1] -= (fwd_rates[42] - rev_rates[42]) * pres_mod[5];
  //sp 4
  sp_rates[4] -= (fwd_rates[42] - rev_rates[42]) * pres_mod[5];
  //sp 5
  sp_rates[5] = (fwd_rates[42] - rev_rates[42]) * pres_mod[5];

  //rxn 43
  //sp 1
  sp_rates[1] -= (fwd_rates[43] - rev_rates[43]);
  //sp 2
  sp_rates[2] += (fwd_rates[43] - rev_rates[43]);
  //sp 5
  sp_rates[5] += (fwd_rates[43] - rev_rates[43]);
  //sp 6
  sp_rates[6] -= (fwd_rates[43] - rev_rates[43]);

  //rxn 44
  //sp 0
  sp_rates[0] += (fwd_rates[44] - rev_rates[44]);
  //sp 1
  sp_rates[1] -= (fwd_rates[44] - rev_rates[44]);
  //sp 3
  sp_rates[3] += (fwd_rates[44] - rev_rates[44]);
  //sp 6
  sp_rates[6] -= (fwd_rates[44] - rev_rates[44]);

  //rxn 45
  //sp 1
  sp_rates[1] -= (fwd_rates[45] - rev_rates[45]);
  //sp 4
  sp_rates[4] += 2.0 * (fwd_rates[45] - rev_rates[45]);
  //sp 6
  sp_rates[6] -= (fwd_rates[45] - rev_rates[45]);

  //rxn 46
  //sp 0
  sp_rates[0] += (fwd_rates[46] - rev_rates[46]);
  //sp 1
  sp_rates[1] -= (fwd_rates[46] - rev_rates[46]);
  //sp 6
  sp_rates[6] += (fwd_rates[46] - rev_rates[46]);
  //sp 7
  sp_rates[7] -= (fwd_rates[46] - rev_rates[46]);

  //rxn 47
  //sp 1
  sp_rates[1] -= (fwd_rates[47] - rev_rates[47]);
  //sp 4
  sp_rates[4] += (fwd_rates[47] - rev_rates[47]);
  //sp 5
  sp_rates[5] += (fwd_rates[47] - rev_rates[47]);
  //sp 7
  sp_rates[7] -= (fwd_rates[47] - rev_rates[47]);

  //rxn 48
  //sp 0
  sp_rates[0] += (fwd_rates[48] - rev_rates[48]);
  //sp 1
  sp_rates[1] -= (fwd_rates[48] - rev_rates[48]);
  //sp 8
  sp_rates[8] = (fwd_rates[48] - rev_rates[48]);
  //sp 9
  sp_rates[9] -= (fwd_rates[48] - rev_rates[48]);

  //rxn 49
  //sp 1
  sp_rates[1] -= (fwd_rates[49] - rev_rates[49]) * pres_mod[6];
  //sp 10
  sp_rates[10] -= (fwd_rates[49] - rev_rates[49]) * pres_mod[6];
  //sp 12
  sp_rates[12] += (fwd_rates[49] - rev_rates[49]) * pres_mod[6];

  //rxn 50
  //sp 0
  sp_rates[0] += (fwd_rates[50] - rev_rates[50]);
  //sp 1
  sp_rates[1] -= (fwd_rates[50] - rev_rates[50]);
  //sp 11
  sp_rates[11] -= (fwd_rates[50] - rev_rates[50]);
  //sp 9
  sp_rates[9] += (fwd_rates[50] - rev_rates[50]);

  //rxn 51
  //sp 1
  sp_rates[1] -= (fwd_rates[51] - rev_rates[51]) * pres_mod[7];
  //sp 12
  sp_rates[12] -= (fwd_rates[51] - rev_rates[51]) * pres_mod[7];
  //sp 13
  sp_rates[13] += (fwd_rates[51] - rev_rates[51]) * pres_mod[7];

  //rxn 52
  //sp 0
  sp_rates[0] += (fwd_rates[52] - rev_rates[52]);
  //sp 1
  sp_rates[1] -= (fwd_rates[52] - rev_rates[52]);
  //sp 12
  sp_rates[12] += (fwd_rates[52] - rev_rates[52]);
  //sp 13
  sp_rates[13] -= (fwd_rates[52] - rev_rates[52]);

  //rxn 53
  //sp 16
  sp_rates[16] -= (fwd_rates[53] - rev_rates[53]) * pres_mod[8];
  //sp 1
  sp_rates[1] -= (fwd_rates[53] - rev_rates[53]) * pres_mod[8];
  //sp 17
  sp_rates[17] += (fwd_rates[53] - rev_rates[53]) * pres_mod[8];

  //rxn 54
  //sp 16
  sp_rates[16] -= (fwd_rates[54] - rev_rates[54]);
  //sp 1
  sp_rates[1] -= (fwd_rates[54] - rev_rates[54]);
  //sp 14
  sp_rates[14] += (fwd_rates[54] - rev_rates[54]);
  //sp 0
  sp_rates[0] += (fwd_rates[54] - rev_rates[54]);

  //rxn 55
  //sp 17
  sp_rates[17] -= (fwd_rates[55] - rev_rates[55]) * pres_mod[9];
  //sp 18
  sp_rates[18] += (fwd_rates[55] - rev_rates[55]) * pres_mod[9];
  //sp 1
  sp_rates[1] -= (fwd_rates[55] - rev_rates[55]) * pres_mod[9];

  //rxn 56
  //sp 17
  sp_rates[17] -= (fwd_rates[56] - rev_rates[56]) * pres_mod[10];
  //sp 19
  sp_rates[19] += (fwd_rates[56] - rev_rates[56]) * pres_mod[10];
  //sp 1
  sp_rates[1] -= (fwd_rates[56] - rev_rates[56]) * pres_mod[10];

  //rxn 57
  //sp 0
  sp_rates[0] += (fwd_rates[57] - rev_rates[57]);
  //sp 17
  sp_rates[17] -= (fwd_rates[57] - rev_rates[57]);
  //sp 16
  sp_rates[16] += (fwd_rates[57] - rev_rates[57]);
  //sp 1
  sp_rates[1] -= (fwd_rates[57] - rev_rates[57]);

  //rxn 58
  //sp 1
  sp_rates[1] -= (fwd_rates[58] - rev_rates[58]) * pres_mod[11];
  //sp 18
  sp_rates[18] -= (fwd_rates[58] - rev_rates[58]) * pres_mod[11];
  //sp 20
  sp_rates[20] += (fwd_rates[58] - rev_rates[58]) * pres_mod[11];

  //rxn 59
  //sp 0
  sp_rates[0] += (fwd_rates[59] - rev_rates[59]);
  //sp 1
  sp_rates[1] -= (fwd_rates[59] - rev_rates[59]);
  //sp 18
  sp_rates[18] -= (fwd_rates[59] - rev_rates[59]);
  //sp 17
  sp_rates[17] += (fwd_rates[59] - rev_rates[59]);

  //rxn 60
  //sp 1
  sp_rates[1] -= (fwd_rates[60] - rev_rates[60]);
  //sp 18
  sp_rates[18] -= (fwd_rates[60] - rev_rates[60]);
  //sp 12
  sp_rates[12] += (fwd_rates[60] - rev_rates[60]);
  //sp 4
  sp_rates[4] += (fwd_rates[60] - rev_rates[60]);

  //rxn 61
  //sp 1
  sp_rates[1] -= (fwd_rates[61] - rev_rates[61]);
  //sp 18
  sp_rates[18] -= (fwd_rates[61] - rev_rates[61]);
  //sp 11
  sp_rates[11] += (fwd_rates[61] - rev_rates[61]);
  //sp 5
  sp_rates[5] += (fwd_rates[61] - rev_rates[61]);

  //rxn 62
  //sp 1
  sp_rates[1] -= (fwd_rates[62] - rev_rates[62]) * pres_mod[12];
  //sp 19
  sp_rates[19] -= (fwd_rates[62] - rev_rates[62]) * pres_mod[12];
  //sp 20
  sp_rates[20] += (fwd_rates[62] - rev_rates[62]) * pres_mod[12];

  //rxn 63
  //sp 18
  sp_rates[18] += (fwd_rates[63] - rev_rates[63]);
  //sp 19
  sp_rates[19] -= (fwd_rates[63] - rev_rates[63]);

  //rxn 64
  //sp 0
  sp_rates[0] += (fwd_rates[64] - rev_rates[64]);
  //sp 1
  sp_rates[1] -= (fwd_rates[64] - rev_rates[64]);
  //sp 19
  sp_rates[19] -= (fwd_rates[64] - rev_rates[64]);
  //sp 17
  sp_rates[17] += (fwd_rates[64] - rev_rates[64]);

  //rxn 65
  //sp 1
  sp_rates[1] -= (fwd_rates[65] - rev_rates[65]);
  //sp 19
  sp_rates[19] -= (fwd_rates[65] - rev_rates[65]);
  //sp 12
  sp_rates[12] += (fwd_rates[65] - rev_rates[65]);
  //sp 4
  sp_rates[4] += (fwd_rates[65] - rev_rates[65]);

  //rxn 66
  //sp 11
  sp_rates[11] += (fwd_rates[66] - rev_rates[66]);
  //sp 1
  sp_rates[1] -= (fwd_rates[66] - rev_rates[66]);
  //sp 19
  sp_rates[19] -= (fwd_rates[66] - rev_rates[66]);
  //sp 5
  sp_rates[5] += (fwd_rates[66] - rev_rates[66]);

  //rxn 67
  //sp 0
  sp_rates[0] += (fwd_rates[67] - rev_rates[67]);
  //sp 1
  sp_rates[1] -= (fwd_rates[67] - rev_rates[67]);
  //sp 18
  sp_rates[18] += (fwd_rates[67] - rev_rates[67]);
  //sp 20
  sp_rates[20] -= (fwd_rates[67] - rev_rates[67]);

  //rxn 68
  //sp 0
  sp_rates[0] += (fwd_rates[68] - rev_rates[68]);
  //sp 1
  sp_rates[1] -= (fwd_rates[68] - rev_rates[68]);
  //sp 19
  sp_rates[19] += (fwd_rates[68] - rev_rates[68]);
  //sp 20
  sp_rates[20] -= (fwd_rates[68] - rev_rates[68]);

  //rxn 69
  //sp 1
  sp_rates[1] -= (fwd_rates[69] - rev_rates[69]) * pres_mod[13];
  //sp 21
  sp_rates[21] -= (fwd_rates[69] - rev_rates[69]) * pres_mod[13];
  //sp 22
  sp_rates[22] += (fwd_rates[69] - rev_rates[69]) * pres_mod[13];

  //rxn 70
  //sp 1
  sp_rates[1] -= (fwd_rates[70] - rev_rates[70]) * pres_mod[14];
  //sp 22
  sp_rates[22] -= (fwd_rates[70] - rev_rates[70]) * pres_mod[14];
  //sp 23
  sp_rates[23] += (fwd_rates[70] - rev_rates[70]) * pres_mod[14];

  //rxn 71
  //sp 24
  sp_rates[24] += (fwd_rates[71] - rev_rates[71]) * pres_mod[15];
  //sp 1
  sp_rates[1] -= (fwd_rates[71] - rev_rates[71]) * pres_mod[15];
  //sp 23
  sp_rates[23] -= (fwd_rates[71] - rev_rates[71]) * pres_mod[15];

  //rxn 72
  //sp 0
  sp_rates[0] += (fwd_rates[72] - rev_rates[72]);
  //sp 1
  sp_rates[1] -= (fwd_rates[72] - rev_rates[72]);
  //sp 22
  sp_rates[22] += (fwd_rates[72] - rev_rates[72]);
  //sp 23
  sp_rates[23] -= (fwd_rates[72] - rev_rates[72]);

  //rxn 73
  //sp 24
  sp_rates[24] -= (fwd_rates[73] - rev_rates[73]) * pres_mod[16];
  //sp 1
  sp_rates[1] -= (fwd_rates[73] - rev_rates[73]) * pres_mod[16];
  //sp 25
  sp_rates[25] += (fwd_rates[73] - rev_rates[73]) * pres_mod[16];

  //rxn 74
  //sp 24
  sp_rates[24] -= (fwd_rates[74] - rev_rates[74]);
  //sp 1
  sp_rates[1] -= (fwd_rates[74] - rev_rates[74]);
  //sp 0
  sp_rates[0] += (fwd_rates[74] - rev_rates[74]);
  //sp 23
  sp_rates[23] += (fwd_rates[74] - rev_rates[74]);

  //rxn 75
  //sp 1
  sp_rates[1] -= (fwd_rates[75] - rev_rates[75]) * pres_mod[17];
  //sp 26
  sp_rates[26] += (fwd_rates[75] - rev_rates[75]) * pres_mod[17];
  //sp 25
  sp_rates[25] -= (fwd_rates[75] - rev_rates[75]) * pres_mod[17];

  //rxn 76
  //sp 0
  sp_rates[0] += (fwd_rates[76] - rev_rates[76]);
  //sp 1
  sp_rates[1] -= (fwd_rates[76] - rev_rates[76]);
  //sp 24
  sp_rates[24] += (fwd_rates[76] - rev_rates[76]);
  //sp 25
  sp_rates[25] -= (fwd_rates[76] - rev_rates[76]);

  //rxn 77
  //sp 0
  sp_rates[0] += (fwd_rates[77] - rev_rates[77]);
  //sp 1
  sp_rates[1] -= (fwd_rates[77] - rev_rates[77]);
  //sp 26
  sp_rates[26] -= (fwd_rates[77] - rev_rates[77]);
  //sp 25
  sp_rates[25] += (fwd_rates[77] - rev_rates[77]);

  //rxn 78
  //sp 11
  sp_rates[11] += (fwd_rates[78] - rev_rates[78]);
  //sp 1
  sp_rates[1] -= (fwd_rates[78] - rev_rates[78]);
  //sp 27
  sp_rates[27] -= (fwd_rates[78] - rev_rates[78]);
  //sp 14
  sp_rates[14] += (fwd_rates[78] - rev_rates[78]);

  //rxn 79
  //sp 0
  sp_rates[0] += (fwd_rates[79] - rev_rates[79]);
  //sp 1
  sp_rates[1] -= (fwd_rates[79] - rev_rates[79]);
  //sp 27
  sp_rates[27] += (fwd_rates[79] - rev_rates[79]);
  //sp 28
  sp_rates[28] -= (fwd_rates[79] - rev_rates[79]);

  //rxn 80
  //sp 1
  sp_rates[1] -= (fwd_rates[80] - rev_rates[80]);
  //sp 28
  sp_rates[28] -= (fwd_rates[80] - rev_rates[80]);
  //sp 14
  sp_rates[14] += (fwd_rates[80] - rev_rates[80]);
  //sp 12
  sp_rates[12] += (fwd_rates[80] - rev_rates[80]);

  //rxn 81
  //sp 28
  sp_rates[28] += (fwd_rates[81] - rev_rates[81]);
  //sp 29
  sp_rates[29] = -(fwd_rates[81] - rev_rates[81]);

  //rxn 82
  //sp 0
  sp_rates[0] -= (fwd_rates[82] - rev_rates[82]) * pres_mod[18];
  //sp 17
  sp_rates[17] += (fwd_rates[82] - rev_rates[82]) * pres_mod[18];
  //sp 14
  sp_rates[14] -= (fwd_rates[82] - rev_rates[82]) * pres_mod[18];

  //rxn 83
  //sp 0
  sp_rates[0] -= (fwd_rates[83] - rev_rates[83]);
  //sp 1
  sp_rates[1] += (fwd_rates[83] - rev_rates[83]);
  //sp 4
  sp_rates[4] -= (fwd_rates[83] - rev_rates[83]);
  //sp 5
  sp_rates[5] += (fwd_rates[83] - rev_rates[83]);

  //rxn 84
  //sp 4
  sp_rates[4] -= 2.0 * (fwd_rates[84] - rev_rates[84]) * pres_mod[19];
  //sp 7
  sp_rates[7] += (fwd_rates[84] - rev_rates[84]) * pres_mod[19];

  //rxn 85
  //sp 2
  sp_rates[2] += (fwd_rates[85] - rev_rates[85]);
  //sp 4
  sp_rates[4] -= 2.0 * (fwd_rates[85] - rev_rates[85]);
  //sp 5
  sp_rates[5] += (fwd_rates[85] - rev_rates[85]);

  //rxn 86
  //sp 3
  sp_rates[3] += (fwd_rates[86] - rev_rates[86]);
  //sp 4
  sp_rates[4] -= (fwd_rates[86] - rev_rates[86]);
  //sp 5
  sp_rates[5] += (fwd_rates[86] - rev_rates[86]);
  //sp 6
  sp_rates[6] -= (fwd_rates[86] - rev_rates[86]);

  //rxn 87
  //sp 4
  sp_rates[4] -= (fwd_rates[87] - rev_rates[87]);
  //sp 5
  sp_rates[5] += (fwd_rates[87] - rev_rates[87]);
  //sp 6
  sp_rates[6] += (fwd_rates[87] - rev_rates[87]);
  //sp 7
  sp_rates[7] -= (fwd_rates[87] - rev_rates[87]);

  //rxn 88
  //sp 4
  sp_rates[4] -= (fwd_rates[88] - rev_rates[88]);
  //sp 5
  sp_rates[5] += (fwd_rates[88] - rev_rates[88]);
  //sp 6
  sp_rates[6] += (fwd_rates[88] - rev_rates[88]);
  //sp 7
  sp_rates[7] -= (fwd_rates[88] - rev_rates[88]);

  //rxn 89
  //sp 8
  sp_rates[8] -= (fwd_rates[89] - rev_rates[89]);
  //sp 1
  sp_rates[1] += (fwd_rates[89] - rev_rates[89]);
  //sp 4
  sp_rates[4] -= (fwd_rates[89] - rev_rates[89]);
  //sp 14
  sp_rates[14] += (fwd_rates[89] - rev_rates[89]);

  //rxn 90
  //sp 16
  sp_rates[16] += (fwd_rates[90] - rev_rates[90]);
  //sp 9
  sp_rates[9] -= (fwd_rates[90] - rev_rates[90]);
  //sp 4
  sp_rates[4] -= (fwd_rates[90] - rev_rates[90]);
  //sp 1
  sp_rates[1] += (fwd_rates[90] - rev_rates[90]);

  //rxn 91
  //sp 17
  sp_rates[17] += (fwd_rates[91] - rev_rates[91]);
  //sp 10
  sp_rates[10] -= (fwd_rates[91] - rev_rates[91]);
  //sp 4
  sp_rates[4] -= (fwd_rates[91] - rev_rates[91]);
  //sp 1
  sp_rates[1] += (fwd_rates[91] - rev_rates[91]);

  //rxn 92
  //sp 9
  sp_rates[9] += (fwd_rates[92] - rev_rates[92]);
  //sp 10
  sp_rates[10] -= (fwd_rates[92] - rev_rates[92]);
  //sp 4
  sp_rates[4] -= (fwd_rates[92] - rev_rates[92]);
  //sp 5
  sp_rates[5] += (fwd_rates[92] - rev_rates[92]);

  //rxn 93
  //sp 17
  sp_rates[17] += (fwd_rates[93] - rev_rates[93]);
  //sp 11
  sp_rates[11] -= (fwd_rates[93] - rev_rates[93]);
  //sp 4
  sp_rates[4] -= (fwd_rates[93] - rev_rates[93]);
  //sp 1
  sp_rates[1] += (fwd_rates[93] - rev_rates[93]);

  //rxn 94
  //sp 4
  sp_rates[4] -= (fwd_rates[94] - rev_rates[94]) * pres_mod[20];
  //sp 12
  sp_rates[12] -= (fwd_rates[94] - rev_rates[94]) * pres_mod[20];
  //sp 20
  sp_rates[20] += (fwd_rates[94] - rev_rates[94]) * pres_mod[20];

  //rxn 95
  //sp 4
  sp_rates[4] -= (fwd_rates[95] - rev_rates[95]);
  //sp 10
  sp_rates[10] += (fwd_rates[95] - rev_rates[95]);
  //sp 12
  sp_rates[12] -= (fwd_rates[95] - rev_rates[95]);
  //sp 5
  sp_rates[5] += (fwd_rates[95] - rev_rates[95]);

  //rxn 96
  //sp 4
  sp_rates[4] -= (fwd_rates[96] - rev_rates[96]);
  //sp 11
  sp_rates[11] += (fwd_rates[96] - rev_rates[96]);
  //sp 12
  sp_rates[12] -= (fwd_rates[96] - rev_rates[96]);
  //sp 5
  sp_rates[5] += (fwd_rates[96] - rev_rates[96]);

  //rxn 97
  //sp 12
  sp_rates[12] += (fwd_rates[97] - rev_rates[97]);
  //sp 4
  sp_rates[4] -= (fwd_rates[97] - rev_rates[97]);
  //sp 13
  sp_rates[13] -= (fwd_rates[97] - rev_rates[97]);
  //sp 5
  sp_rates[5] += (fwd_rates[97] - rev_rates[97]);

  //rxn 98
  //sp 1
  sp_rates[1] += (fwd_rates[98] - rev_rates[98]);
  //sp 4
  sp_rates[4] -= (fwd_rates[98] - rev_rates[98]);
  //sp 14
  sp_rates[14] -= (fwd_rates[98] - rev_rates[98]);
  //sp 15
  sp_rates[15] += (fwd_rates[98] - rev_rates[98]);

  //rxn 99
  //sp 16
  sp_rates[16] -= (fwd_rates[99] - rev_rates[99]);
  //sp 4
  sp_rates[4] -= (fwd_rates[99] - rev_rates[99]);
  //sp 5
  sp_rates[5] += (fwd_rates[99] - rev_rates[99]);
  //sp 14
  sp_rates[14] += (fwd_rates[99] - rev_rates[99]);

  //rxn 100
  //sp 16
  sp_rates[16] += (fwd_rates[100] - rev_rates[100]);
  //sp 17
  sp_rates[17] -= (fwd_rates[100] - rev_rates[100]);
  //sp 4
  sp_rates[4] -= (fwd_rates[100] - rev_rates[100]);
  //sp 5
  sp_rates[5] += (fwd_rates[100] - rev_rates[100]);

  //rxn 101
  //sp 17
  sp_rates[17] += (fwd_rates[101] - rev_rates[101]);
  //sp 18
  sp_rates[18] -= (fwd_rates[101] - rev_rates[101]);
  //sp 4
  sp_rates[4] -= (fwd_rates[101] - rev_rates[101]);
  //sp 5
  sp_rates[5] += (fwd_rates[101] - rev_rates[101]);

  //rxn 102
  //sp 17
  sp_rates[17] += (fwd_rates[102] - rev_rates[102]);
  //sp 19
  sp_rates[19] -= (fwd_rates[102] - rev_rates[102]);
  //sp 4
  sp_rates[4] -= (fwd_rates[102] - rev_rates[102]);
  //sp 5
  sp_rates[5] += (fwd_rates[102] - rev_rates[102]);

  //rxn 103
  //sp 4
  sp_rates[4] -= (fwd_rates[103] - rev_rates[103]);
  //sp 18
  sp_rates[18] += (fwd_rates[103] - rev_rates[103]);
  //sp 20
  sp_rates[20] -= (fwd_rates[103] - rev_rates[103]);
  //sp 5
  sp_rates[5] += (fwd_rates[103] - rev_rates[103]);

  //rxn 104
  //sp 4
  sp_rates[4] -= (fwd_rates[104] - rev_rates[104]);
  //sp 19
  sp_rates[19] += (fwd_rates[104] - rev_rates[104]);
  //sp 20
  sp_rates[20] -= (fwd_rates[104] - rev_rates[104]);
  //sp 5
  sp_rates[5] += (fwd_rates[104] - rev_rates[104]);

  //rxn 105
  //sp 1
  sp_rates[1] += (fwd_rates[105] - rev_rates[105]);
  //sp 27
  sp_rates[27] += (fwd_rates[105] - rev_rates[105]);
  //sp 4
  sp_rates[4] -= (fwd_rates[105] - rev_rates[105]);
  //sp 21
  sp_rates[21] -= (fwd_rates[105] - rev_rates[105]);

  //rxn 106
  //sp 1
  sp_rates[1] += (fwd_rates[106] - rev_rates[106]);
  //sp 4
  sp_rates[4] -= (fwd_rates[106] - rev_rates[106]);
  //sp 22
  sp_rates[22] -= (fwd_rates[106] - rev_rates[106]);
  //sp 28
  sp_rates[28] += (fwd_rates[106] - rev_rates[106]);

  //rxn 107
  //sp 1
  sp_rates[1] += (fwd_rates[107] - rev_rates[107]);
  //sp 4
  sp_rates[4] -= (fwd_rates[107] - rev_rates[107]);
  //sp 29
  sp_rates[29] += (fwd_rates[107] - rev_rates[107]);
  //sp 22
  sp_rates[22] -= (fwd_rates[107] - rev_rates[107]);

  //rxn 108
  //sp 4
  sp_rates[4] -= (fwd_rates[108] - rev_rates[108]);
  //sp 21
  sp_rates[21] += (fwd_rates[108] - rev_rates[108]);
  //sp 22
  sp_rates[22] -= (fwd_rates[108] - rev_rates[108]);
  //sp 5
  sp_rates[5] += (fwd_rates[108] - rev_rates[108]);

  //rxn 109
  //sp 12
  sp_rates[12] += (fwd_rates[109] - rev_rates[109]);
  //sp 4
  sp_rates[4] -= (fwd_rates[109] - rev_rates[109]);
  //sp 14
  sp_rates[14] += (fwd_rates[109] - rev_rates[109]);
  //sp 22
  sp_rates[22] -= (fwd_rates[109] - rev_rates[109]);

  //rxn 110
  //sp 4
  sp_rates[4] -= (fwd_rates[110] - rev_rates[110]);
  //sp 5
  sp_rates[5] += (fwd_rates[110] - rev_rates[110]);
  //sp 22
  sp_rates[22] += (fwd_rates[110] - rev_rates[110]);
  //sp 23
  sp_rates[23] -= (fwd_rates[110] - rev_rates[110]);

  //rxn 111
  //sp 24
  sp_rates[24] -= (fwd_rates[111] - rev_rates[111]);
  //sp 4
  sp_rates[4] -= (fwd_rates[111] - rev_rates[111]);
  //sp 5
  sp_rates[5] += (fwd_rates[111] - rev_rates[111]);
  //sp 23
  sp_rates[23] += (fwd_rates[111] - rev_rates[111]);

  //rxn 112
  //sp 25
  sp_rates[25] += (fwd_rates[112] - rev_rates[112]);
  //sp 26
  sp_rates[26] -= (fwd_rates[112] - rev_rates[112]);
  //sp 4
  sp_rates[4] -= (fwd_rates[112] - rev_rates[112]);
  //sp 5
  sp_rates[5] += (fwd_rates[112] - rev_rates[112]);

  //rxn 113
  //sp 4
  sp_rates[4] -= (fwd_rates[113] - rev_rates[113]);
  //sp 27
  sp_rates[27] += (fwd_rates[113] - rev_rates[113]);
  //sp 28
  sp_rates[28] -= (fwd_rates[113] - rev_rates[113]);
  //sp 5
  sp_rates[5] += (fwd_rates[113] - rev_rates[113]);

  //rxn 114
  //sp 3
  sp_rates[3] += (fwd_rates[114] - rev_rates[114]);
  //sp 6
  sp_rates[6] -= 2.0 * (fwd_rates[114] - rev_rates[114]);
  //sp 7
  sp_rates[7] += (fwd_rates[114] - rev_rates[114]);

  //rxn 115
  //sp 3
  sp_rates[3] += (fwd_rates[115] - rev_rates[115]);
  //sp 6
  sp_rates[6] -= 2.0 * (fwd_rates[115] - rev_rates[115]);
  //sp 7
  sp_rates[7] += (fwd_rates[115] - rev_rates[115]);

  //rxn 116
  //sp 17
  sp_rates[17] += (fwd_rates[116] - rev_rates[116]);
  //sp 10
  sp_rates[10] -= (fwd_rates[116] - rev_rates[116]);
  //sp 4
  sp_rates[4] += (fwd_rates[116] - rev_rates[116]);
  //sp 6
  sp_rates[6] -= (fwd_rates[116] - rev_rates[116]);

  //rxn 117
  //sp 3
  sp_rates[3] += (fwd_rates[117] - rev_rates[117]);
  //sp 12
  sp_rates[12] -= (fwd_rates[117] - rev_rates[117]);
  //sp 13
  sp_rates[13] += (fwd_rates[117] - rev_rates[117]);
  //sp 6
  sp_rates[6] -= (fwd_rates[117] - rev_rates[117]);

  //rxn 118
  //sp 4
  sp_rates[4] += (fwd_rates[118] - rev_rates[118]);
  //sp 19
  sp_rates[19] += (fwd_rates[118] - rev_rates[118]);
  //sp 12
  sp_rates[12] -= (fwd_rates[118] - rev_rates[118]);
  //sp 6
  sp_rates[6] -= (fwd_rates[118] - rev_rates[118]);

  //rxn 119
  //sp 4
  sp_rates[4] += (fwd_rates[119] - rev_rates[119]);
  //sp 6
  sp_rates[6] -= (fwd_rates[119] - rev_rates[119]);
  //sp 14
  sp_rates[14] -= (fwd_rates[119] - rev_rates[119]);
  //sp 15
  sp_rates[15] += (fwd_rates[119] - rev_rates[119]);

  //rxn 120
  //sp 16
  sp_rates[16] += (fwd_rates[120] - rev_rates[120]);
  //sp 17
  sp_rates[17] -= (fwd_rates[120] - rev_rates[120]);
  //sp 6
  sp_rates[6] -= (fwd_rates[120] - rev_rates[120]);
  //sp 7
  sp_rates[7] += (fwd_rates[120] - rev_rates[120]);

  //rxn 121
  //sp 8
  sp_rates[8] -= (fwd_rates[121] - rev_rates[121]);
  //sp 2
  sp_rates[2] += (fwd_rates[121] - rev_rates[121]);
  //sp 3
  sp_rates[3] -= (fwd_rates[121] - rev_rates[121]);
  //sp 14
  sp_rates[14] += (fwd_rates[121] - rev_rates[121]);

  //rxn 122
  //sp 8
  sp_rates[8] -= (fwd_rates[122] - rev_rates[122]);
  //sp 1
  sp_rates[1] += (fwd_rates[122] - rev_rates[122]);
  //sp 10
  sp_rates[10] -= (fwd_rates[122] - rev_rates[122]);
  //sp 21
  sp_rates[21] += (fwd_rates[122] - rev_rates[122]);

  //rxn 123
  //sp 8
  sp_rates[8] -= (fwd_rates[123] - rev_rates[123]);
  //sp 1
  sp_rates[1] += (fwd_rates[123] - rev_rates[123]);
  //sp 12
  sp_rates[12] -= (fwd_rates[123] - rev_rates[123]);
  //sp 22
  sp_rates[22] += (fwd_rates[123] - rev_rates[123]);

  //rxn 124
  //sp 16
  sp_rates[16] += (fwd_rates[124] - rev_rates[124]);
  //sp 9
  sp_rates[9] -= (fwd_rates[124] - rev_rates[124]);
  //sp 2
  sp_rates[2] += (fwd_rates[124] - rev_rates[124]);
  //sp 3
  sp_rates[3] -= (fwd_rates[124] - rev_rates[124]);

  //rxn 125
  //sp 0
  sp_rates[0] -= (fwd_rates[125] - rev_rates[125]);
  //sp 9
  sp_rates[9] -= (fwd_rates[125] - rev_rates[125]);
  //sp 10
  sp_rates[10] += (fwd_rates[125] - rev_rates[125]);
  //sp 1
  sp_rates[1] += (fwd_rates[125] - rev_rates[125]);

  //rxn 126
  //sp 9
  sp_rates[9] -= (fwd_rates[126] - rev_rates[126]);
  //sp 1
  sp_rates[1] += (fwd_rates[126] - rev_rates[126]);
  //sp 5
  sp_rates[5] -= (fwd_rates[126] - rev_rates[126]);
  //sp 17
  sp_rates[17] += (fwd_rates[126] - rev_rates[126]);

  //rxn 127
  //sp 9
  sp_rates[9] -= (fwd_rates[127] - rev_rates[127]);
  //sp 10
  sp_rates[10] -= (fwd_rates[127] - rev_rates[127]);
  //sp 22
  sp_rates[22] += (fwd_rates[127] - rev_rates[127]);
  //sp 1
  sp_rates[1] += (fwd_rates[127] - rev_rates[127]);

  //rxn 128
  //sp 9
  sp_rates[9] -= (fwd_rates[128] - rev_rates[128]);
  //sp 23
  sp_rates[23] += (fwd_rates[128] - rev_rates[128]);
  //sp 12
  sp_rates[12] -= (fwd_rates[128] - rev_rates[128]);
  //sp 1
  sp_rates[1] += (fwd_rates[128] - rev_rates[128]);

  //rxn 129
  //sp 24
  sp_rates[24] += (fwd_rates[129] - rev_rates[129]);
  //sp 9
  sp_rates[9] -= (fwd_rates[129] - rev_rates[129]);
  //sp 13
  sp_rates[13] -= (fwd_rates[129] - rev_rates[129]);
  //sp 1
  sp_rates[1] += (fwd_rates[129] - rev_rates[129]);

  //rxn 130
  //sp 9
  sp_rates[9] -= (fwd_rates[130] - rev_rates[130]) * pres_mod[21];
  //sp 27
  sp_rates[27] += (fwd_rates[130] - rev_rates[130]) * pres_mod[21];
  //sp 14
  sp_rates[14] -= (fwd_rates[130] - rev_rates[130]) * pres_mod[21];

  //rxn 131
  //sp 16
  sp_rates[16] += (fwd_rates[131] - rev_rates[131]);
  //sp 9
  sp_rates[9] -= (fwd_rates[131] - rev_rates[131]);
  //sp 14
  sp_rates[14] += (fwd_rates[131] - rev_rates[131]);
  //sp 15
  sp_rates[15] -= (fwd_rates[131] - rev_rates[131]);

  //rxn 132
  //sp 17
  sp_rates[17] -= (fwd_rates[132] - rev_rates[132]);
  //sp 1
  sp_rates[1] += (fwd_rates[132] - rev_rates[132]);
  //sp 28
  sp_rates[28] += (fwd_rates[132] - rev_rates[132]);
  //sp 9
  sp_rates[9] -= (fwd_rates[132] - rev_rates[132]);

  //rxn 133
  //sp 9
  sp_rates[9] -= (fwd_rates[133] - rev_rates[133]);
  //sp 27
  sp_rates[27] -= (fwd_rates[133] - rev_rates[133]);
  //sp 22
  sp_rates[22] += (fwd_rates[133] - rev_rates[133]);
  //sp 14
  sp_rates[14] += (fwd_rates[133] - rev_rates[133]);

  //rxn 134
  //sp 1
  sp_rates[1] += fwd_rates[134];
  //sp 10
  sp_rates[10] -= fwd_rates[134];
  //sp 3
  sp_rates[3] -= fwd_rates[134];
  //sp 4
  sp_rates[4] += fwd_rates[134];
  //sp 14
  sp_rates[14] += fwd_rates[134];

  //rxn 135
  //sp 0
  sp_rates[0] -= (fwd_rates[135] - rev_rates[134]);
  //sp 1
  sp_rates[1] += (fwd_rates[135] - rev_rates[134]);
  //sp 10
  sp_rates[10] -= (fwd_rates[135] - rev_rates[134]);
  //sp 12
  sp_rates[12] += (fwd_rates[135] - rev_rates[134]);

  //rxn 136
  //sp 0
  sp_rates[0] += (fwd_rates[136] - rev_rates[135]);
  //sp 10
  sp_rates[10] -= 2.0 * (fwd_rates[136] - rev_rates[135]);
  //sp 22
  sp_rates[22] += (fwd_rates[136] - rev_rates[135]);

  //rxn 137
  //sp 24
  sp_rates[24] += (fwd_rates[137] - rev_rates[136]);
  //sp 1
  sp_rates[1] += (fwd_rates[137] - rev_rates[136]);
  //sp 10
  sp_rates[10] -= (fwd_rates[137] - rev_rates[136]);
  //sp 12
  sp_rates[12] -= (fwd_rates[137] - rev_rates[136]);

  //rxn 138
  //sp 10
  sp_rates[10] -= (fwd_rates[138] - rev_rates[137]);
  //sp 12
  sp_rates[12] += 2.0 * (fwd_rates[138] - rev_rates[137]);
  //sp 13
  sp_rates[13] -= (fwd_rates[138] - rev_rates[137]);

  //rxn 139
  //sp 10
  sp_rates[10] -= (fwd_rates[139] - rev_rates[138]) * pres_mod[22];
  //sp 28
  sp_rates[28] += (fwd_rates[139] - rev_rates[138]) * pres_mod[22];
  //sp 14
  sp_rates[14] -= (fwd_rates[139] - rev_rates[138]) * pres_mod[22];

  //rxn 140
  //sp 10
  sp_rates[10] -= (fwd_rates[140] - rev_rates[139]);
  //sp 27
  sp_rates[27] -= (fwd_rates[140] - rev_rates[139]);
  //sp 14
  sp_rates[14] += (fwd_rates[140] - rev_rates[139]);
  //sp 23
  sp_rates[23] += (fwd_rates[140] - rev_rates[139]);

  //rxn 141
  //sp 10
  sp_rates[10] += (fwd_rates[141] - rev_rates[140]);
  //sp 11
  sp_rates[11] -= (fwd_rates[141] - rev_rates[140]);

  //rxn 142
  //sp 10
  sp_rates[10] += (fwd_rates[142] - rev_rates[141]);
  //sp 11
  sp_rates[11] -= (fwd_rates[142] - rev_rates[141]);

  //rxn 143
  //sp 3
  sp_rates[3] -= (fwd_rates[143] - rev_rates[142]);
  //sp 1
  sp_rates[1] += (fwd_rates[143] - rev_rates[142]);
  //sp 11
  sp_rates[11] -= (fwd_rates[143] - rev_rates[142]);
  //sp 4
  sp_rates[4] += (fwd_rates[143] - rev_rates[142]);
  //sp 14
  sp_rates[14] += (fwd_rates[143] - rev_rates[142]);

  //rxn 144
  //sp 3
  sp_rates[3] -= (fwd_rates[144] - rev_rates[143]);
  //sp 11
  sp_rates[11] -= (fwd_rates[144] - rev_rates[143]);
  //sp 5
  sp_rates[5] += (fwd_rates[144] - rev_rates[143]);
  //sp 14
  sp_rates[14] += (fwd_rates[144] - rev_rates[143]);

  //rxn 145
  //sp 0
  sp_rates[0] -= (fwd_rates[145] - rev_rates[144]);
  //sp 1
  sp_rates[1] += (fwd_rates[145] - rev_rates[144]);
  //sp 11
  sp_rates[11] -= (fwd_rates[145] - rev_rates[144]);
  //sp 12
  sp_rates[12] += (fwd_rates[145] - rev_rates[144]);

  //rxn 146
  //sp 11
  sp_rates[11] -= (fwd_rates[146] - rev_rates[145]) * pres_mod[23];
  //sp 20
  sp_rates[20] += (fwd_rates[146] - rev_rates[145]) * pres_mod[23];
  //sp 5
  sp_rates[5] -= (fwd_rates[146] - rev_rates[145]) * pres_mod[23];

  //rxn 147
  //sp 10
  sp_rates[10] += (fwd_rates[147] - rev_rates[146]);
  //sp 11
  sp_rates[11] -= (fwd_rates[147] - rev_rates[146]);

  //rxn 148
  //sp 24
  sp_rates[24] += (fwd_rates[148] - rev_rates[147]);
  //sp 1
  sp_rates[1] += (fwd_rates[148] - rev_rates[147]);
  //sp 11
  sp_rates[11] -= (fwd_rates[148] - rev_rates[147]);
  //sp 12
  sp_rates[12] -= (fwd_rates[148] - rev_rates[147]);

  //rxn 149
  //sp 11
  sp_rates[11] -= (fwd_rates[149] - rev_rates[148]);
  //sp 12
  sp_rates[12] += 2.0 * (fwd_rates[149] - rev_rates[148]);
  //sp 13
  sp_rates[13] -= (fwd_rates[149] - rev_rates[148]);

  //rxn 150
  //sp 10
  sp_rates[10] += (fwd_rates[150] - rev_rates[149]);
  //sp 11
  sp_rates[11] -= (fwd_rates[150] - rev_rates[149]);

  //rxn 151
  //sp 10
  sp_rates[10] += (fwd_rates[151] - rev_rates[150]);
  //sp 11
  sp_rates[11] -= (fwd_rates[151] - rev_rates[150]);

  //rxn 152
  //sp 17
  sp_rates[17] += (fwd_rates[152] - rev_rates[151]);
  //sp 11
  sp_rates[11] -= (fwd_rates[152] - rev_rates[151]);
  //sp 14
  sp_rates[14] += (fwd_rates[152] - rev_rates[151]);
  //sp 15
  sp_rates[15] -= (fwd_rates[152] - rev_rates[151]);

  //rxn 153
  //sp 25
  sp_rates[25] += (fwd_rates[153] - rev_rates[152]);
  //sp 26
  sp_rates[26] -= (fwd_rates[153] - rev_rates[152]);
  //sp 11
  sp_rates[11] -= (fwd_rates[153] - rev_rates[152]);
  //sp 12
  sp_rates[12] += (fwd_rates[153] - rev_rates[152]);

  //rxn 154
  //sp 19
  sp_rates[19] += (fwd_rates[154] - rev_rates[153]);
  //sp 2
  sp_rates[2] += (fwd_rates[154] - rev_rates[153]);
  //sp 3
  sp_rates[3] -= (fwd_rates[154] - rev_rates[153]);
  //sp 12
  sp_rates[12] -= (fwd_rates[154] - rev_rates[153]);

  //rxn 155
  //sp 17
  sp_rates[17] += (fwd_rates[155] - rev_rates[154]);
  //sp 3
  sp_rates[3] -= (fwd_rates[155] - rev_rates[154]);
  //sp 12
  sp_rates[12] -= (fwd_rates[155] - rev_rates[154]);
  //sp 4
  sp_rates[4] += (fwd_rates[155] - rev_rates[154]);

  //rxn 156
  //sp 12
  sp_rates[12] -= (fwd_rates[156] - rev_rates[155]);
  //sp 13
  sp_rates[13] += (fwd_rates[156] - rev_rates[155]);
  //sp 6
  sp_rates[6] += (fwd_rates[156] - rev_rates[155]);
  //sp 7
  sp_rates[7] -= (fwd_rates[156] - rev_rates[155]);

  //rxn 157
  //sp 26
  sp_rates[26] += (fwd_rates[157] - rev_rates[156]) * pres_mod[24];
  //sp 12
  sp_rates[12] -= 2.0 * (fwd_rates[157] - rev_rates[156]) * pres_mod[24];

  //rxn 158
  //sp 1
  sp_rates[1] += (fwd_rates[158] - rev_rates[157]);
  //sp 12
  sp_rates[12] -= 2.0 * (fwd_rates[158] - rev_rates[157]);
  //sp 25
  sp_rates[25] += (fwd_rates[158] - rev_rates[157]);

  //rxn 159
  //sp 16
  sp_rates[16] -= (fwd_rates[159] - rev_rates[158]);
  //sp 12
  sp_rates[12] -= (fwd_rates[159] - rev_rates[158]);
  //sp 13
  sp_rates[13] += (fwd_rates[159] - rev_rates[158]);
  //sp 14
  sp_rates[14] += (fwd_rates[159] - rev_rates[158]);

  //rxn 160
  //sp 16
  sp_rates[16] += (fwd_rates[160] - rev_rates[159]);
  //sp 17
  sp_rates[17] -= (fwd_rates[160] - rev_rates[159]);
  //sp 12
  sp_rates[12] -= (fwd_rates[160] - rev_rates[159]);
  //sp 13
  sp_rates[13] += (fwd_rates[160] - rev_rates[159]);

  //rxn 161
  //sp 12
  sp_rates[12] -= (fwd_rates[161] - rev_rates[160]);
  //sp 18
  sp_rates[18] += (fwd_rates[161] - rev_rates[160]);
  //sp 20
  sp_rates[20] -= (fwd_rates[161] - rev_rates[160]);
  //sp 13
  sp_rates[13] += (fwd_rates[161] - rev_rates[160]);

  //rxn 162
  //sp 12
  sp_rates[12] -= (fwd_rates[162] - rev_rates[161]);
  //sp 19
  sp_rates[19] += (fwd_rates[162] - rev_rates[161]);
  //sp 20
  sp_rates[20] -= (fwd_rates[162] - rev_rates[161]);
  //sp 13
  sp_rates[13] += (fwd_rates[162] - rev_rates[161]);

  //rxn 163
  //sp 24
  sp_rates[24] -= (fwd_rates[163] - rev_rates[162]);
  //sp 12
  sp_rates[12] -= (fwd_rates[163] - rev_rates[162]);
  //sp 13
  sp_rates[13] += (fwd_rates[163] - rev_rates[162]);
  //sp 23
  sp_rates[23] += (fwd_rates[163] - rev_rates[162]);

  //rxn 164
  //sp 25
  sp_rates[25] += (fwd_rates[164] - rev_rates[163]);
  //sp 26
  sp_rates[26] -= (fwd_rates[164] - rev_rates[163]);
  //sp 12
  sp_rates[12] -= (fwd_rates[164] - rev_rates[163]);
  //sp 13
  sp_rates[13] += (fwd_rates[164] - rev_rates[163]);

  //rxn 165
  //sp 16
  sp_rates[16] -= (fwd_rates[165] - rev_rates[164]);
  //sp 1
  sp_rates[1] += (fwd_rates[165] - rev_rates[164]);
  //sp 14
  sp_rates[14] += (fwd_rates[165] - rev_rates[164]);

  //rxn 166
  //sp 16
  sp_rates[16] -= (fwd_rates[166] - rev_rates[165]) * pres_mod[25];
  //sp 1
  sp_rates[1] += (fwd_rates[166] - rev_rates[165]) * pres_mod[25];
  //sp 14
  sp_rates[14] += (fwd_rates[166] - rev_rates[165]) * pres_mod[25];

  //rxn 167
  //sp 16
  sp_rates[16] -= (fwd_rates[167] - rev_rates[166]);
  //sp 3
  sp_rates[3] -= (fwd_rates[167] - rev_rates[166]);
  //sp 6
  sp_rates[6] += (fwd_rates[167] - rev_rates[166]);
  //sp 14
  sp_rates[14] += (fwd_rates[167] - rev_rates[166]);

  //rxn 168
  //sp 17
  sp_rates[17] += (fwd_rates[168] - rev_rates[167]);
  //sp 18
  sp_rates[18] -= (fwd_rates[168] - rev_rates[167]);
  //sp 3
  sp_rates[3] -= (fwd_rates[168] - rev_rates[167]);
  //sp 6
  sp_rates[6] += (fwd_rates[168] - rev_rates[167]);

  //rxn 169
  //sp 3
  sp_rates[3] -= (fwd_rates[169] - rev_rates[168]);
  //sp 17
  sp_rates[17] += (fwd_rates[169] - rev_rates[168]);
  //sp 19
  sp_rates[19] -= (fwd_rates[169] - rev_rates[168]);
  //sp 6
  sp_rates[6] += (fwd_rates[169] - rev_rates[168]);

  //rxn 170
  //sp 16
  sp_rates[16] += (fwd_rates[170] - rev_rates[169]);
  //sp 3
  sp_rates[3] -= (fwd_rates[170] - rev_rates[169]);
  //sp 21
  sp_rates[21] -= (fwd_rates[170] - rev_rates[169]);
  //sp 14
  sp_rates[14] += (fwd_rates[170] - rev_rates[169]);

  //rxn 171
  //sp 0
  sp_rates[0] -= (fwd_rates[171] - rev_rates[170]);
  //sp 1
  sp_rates[1] += (fwd_rates[171] - rev_rates[170]);
  //sp 21
  sp_rates[21] -= (fwd_rates[171] - rev_rates[170]);
  //sp 22
  sp_rates[22] += (fwd_rates[171] - rev_rates[170]);

  //rxn 172
  //sp 16
  sp_rates[16] += (fwd_rates[172] - rev_rates[171]);
  //sp 17
  sp_rates[17] += (fwd_rates[172] - rev_rates[171]);
  //sp 3
  sp_rates[3] -= (fwd_rates[172] - rev_rates[171]);
  //sp 23
  sp_rates[23] -= (fwd_rates[172] - rev_rates[171]);

  //rxn 173
  //sp 24
  sp_rates[24] -= (fwd_rates[173] - rev_rates[172]) * pres_mod[26];
  //sp 0
  sp_rates[0] += (fwd_rates[173] - rev_rates[172]) * pres_mod[26];
  //sp 22
  sp_rates[22] += (fwd_rates[173] - rev_rates[172]) * pres_mod[26];

  //rxn 174
  //sp 24
  sp_rates[24] += (fwd_rates[174] - rev_rates[173]);
  //sp 25
  sp_rates[25] -= (fwd_rates[174] - rev_rates[173]);
  //sp 3
  sp_rates[3] -= (fwd_rates[174] - rev_rates[173]);
  //sp 6
  sp_rates[6] += (fwd_rates[174] - rev_rates[173]);

  //rxn 175
  //sp 3
  sp_rates[3] -= (fwd_rates[175] - rev_rates[174]);
  //sp 27
  sp_rates[27] -= (fwd_rates[175] - rev_rates[174]);
  //sp 4
  sp_rates[4] += (fwd_rates[175] - rev_rates[174]);
  //sp 14
  sp_rates[14] += 2.0 * (fwd_rates[175] - rev_rates[174]);

  //rxn 176
  //sp 27
  sp_rates[27] -= 2.0 * (fwd_rates[176] - rev_rates[175]);
  //sp 22
  sp_rates[22] += (fwd_rates[176] - rev_rates[175]);
  //sp 14
  sp_rates[14] += 2.0 * (fwd_rates[176] - rev_rates[175]);

  //rxn 177
  //sp 2
  sp_rates[2] += (fwd_rates[177] - rev_rates[176]);
  //sp 35
  sp_rates[35] = -(fwd_rates[177] - rev_rates[176]);
  //sp 47
  (*dy_N) = (fwd_rates[177] - rev_rates[176]);
  //sp 30
  sp_rates[30] = -(fwd_rates[177] - rev_rates[176]);

  //rxn 178
  //sp 35
  sp_rates[35] += (fwd_rates[178] - rev_rates[177]);
  //sp 2
  sp_rates[2] += (fwd_rates[178] - rev_rates[177]);
  //sp 3
  sp_rates[3] -= (fwd_rates[178] - rev_rates[177]);
  //sp 30
  sp_rates[30] -= (fwd_rates[178] - rev_rates[177]);

  //rxn 179
  //sp 1
  sp_rates[1] += (fwd_rates[179] - rev_rates[178]);
  //sp 35
  sp_rates[35] += (fwd_rates[179] - rev_rates[178]);
  //sp 4
  sp_rates[4] -= (fwd_rates[179] - rev_rates[178]);
  //sp 30
  sp_rates[30] -= (fwd_rates[179] - rev_rates[178]);

  //rxn 180
  //sp 2
  sp_rates[2] -= (fwd_rates[180] - rev_rates[179]);
  //sp 3
  sp_rates[3] += (fwd_rates[180] - rev_rates[179]);
  //sp 47
  (*dy_N) += (fwd_rates[180] - rev_rates[179]);
  //sp 37
  sp_rates[37] = -(fwd_rates[180] - rev_rates[179]);

  //rxn 181
  //sp 2
  sp_rates[2] -= (fwd_rates[181] - rev_rates[180]);
  //sp 35
  sp_rates[35] += 2.0 * (fwd_rates[181] - rev_rates[180]);
  //sp 37
  sp_rates[37] -= (fwd_rates[181] - rev_rates[180]);

  //rxn 182
  //sp 1
  sp_rates[1] -= (fwd_rates[182] - rev_rates[181]);
  //sp 47
  (*dy_N) += (fwd_rates[182] - rev_rates[181]);
  //sp 37
  sp_rates[37] -= (fwd_rates[182] - rev_rates[181]);
  //sp 4
  sp_rates[4] += (fwd_rates[182] - rev_rates[181]);

  //rxn 183
  //sp 47
  (*dy_N) += (fwd_rates[183] - rev_rates[182]);
  //sp 4
  sp_rates[4] -= (fwd_rates[183] - rev_rates[182]);
  //sp 37
  sp_rates[37] -= (fwd_rates[183] - rev_rates[182]);
  //sp 6
  sp_rates[6] += (fwd_rates[183] - rev_rates[182]);

  //rxn 184
  //sp 2
  sp_rates[2] += (fwd_rates[184] - rev_rates[183]) * pres_mod[27];
  //sp 47
  (*dy_N) += (fwd_rates[184] - rev_rates[183]) * pres_mod[27];
  //sp 37
  sp_rates[37] -= (fwd_rates[184] - rev_rates[183]) * pres_mod[27];

  //rxn 185
  //sp 4
  sp_rates[4] += (fwd_rates[185] - rev_rates[184]);
  //sp 35
  sp_rates[35] -= (fwd_rates[185] - rev_rates[184]);
  //sp 36
  sp_rates[36] = (fwd_rates[185] - rev_rates[184]);
  //sp 6
  sp_rates[6] -= (fwd_rates[185] - rev_rates[184]);

  //rxn 186
  //sp 2
  sp_rates[2] -= (fwd_rates[186] - rev_rates[185]) * pres_mod[28];
  //sp 35
  sp_rates[35] -= (fwd_rates[186] - rev_rates[185]) * pres_mod[28];
  //sp 36
  sp_rates[36] += (fwd_rates[186] - rev_rates[185]) * pres_mod[28];

  //rxn 187
  //sp 35
  sp_rates[35] += (fwd_rates[187] - rev_rates[186]);
  //sp 2
  sp_rates[2] -= (fwd_rates[187] - rev_rates[186]);
  //sp 3
  sp_rates[3] += (fwd_rates[187] - rev_rates[186]);
  //sp 36
  sp_rates[36] -= (fwd_rates[187] - rev_rates[186]);

  //rxn 188
  //sp 1
  sp_rates[1] -= (fwd_rates[188] - rev_rates[187]);
  //sp 35
  sp_rates[35] += (fwd_rates[188] - rev_rates[187]);
  //sp 36
  sp_rates[36] -= (fwd_rates[188] - rev_rates[187]);
  //sp 4
  sp_rates[4] += (fwd_rates[188] - rev_rates[187]);

  //rxn 189
  //sp 1
  sp_rates[1] += (fwd_rates[189] - rev_rates[188]);
  //sp 2
  sp_rates[2] -= (fwd_rates[189] - rev_rates[188]);
  //sp 35
  sp_rates[35] += (fwd_rates[189] - rev_rates[188]);
  //sp 31
  sp_rates[31] = -(fwd_rates[189] - rev_rates[188]);

  //rxn 190
  //sp 0
  sp_rates[0] += (fwd_rates[190] - rev_rates[189]);
  //sp 1
  sp_rates[1] -= (fwd_rates[190] - rev_rates[189]);
  //sp 30
  sp_rates[30] += (fwd_rates[190] - rev_rates[189]);
  //sp 31
  sp_rates[31] -= (fwd_rates[190] - rev_rates[189]);

  //rxn 191
  //sp 1
  sp_rates[1] += (fwd_rates[191] - rev_rates[190]);
  //sp 4
  sp_rates[4] -= (fwd_rates[191] - rev_rates[190]);
  //sp 38
  sp_rates[38] = (fwd_rates[191] - rev_rates[190]);
  //sp 31
  sp_rates[31] -= (fwd_rates[191] - rev_rates[190]);

  //rxn 192
  //sp 4
  sp_rates[4] -= (fwd_rates[192] - rev_rates[191]);
  //sp 5
  sp_rates[5] += (fwd_rates[192] - rev_rates[191]);
  //sp 30
  sp_rates[30] += (fwd_rates[192] - rev_rates[191]);
  //sp 31
  sp_rates[31] -= (fwd_rates[192] - rev_rates[191]);

  //rxn 193
  //sp 2
  sp_rates[2] += (fwd_rates[193] - rev_rates[192]);
  //sp 3
  sp_rates[3] -= (fwd_rates[193] - rev_rates[192]);
  //sp 38
  sp_rates[38] += (fwd_rates[193] - rev_rates[192]);
  //sp 31
  sp_rates[31] -= (fwd_rates[193] - rev_rates[192]);

  //rxn 194
  //sp 35
  sp_rates[35] += (fwd_rates[194] - rev_rates[193]);
  //sp 3
  sp_rates[3] -= (fwd_rates[194] - rev_rates[193]);
  //sp 4
  sp_rates[4] += (fwd_rates[194] - rev_rates[193]);
  //sp 31
  sp_rates[31] -= (fwd_rates[194] - rev_rates[193]);

  //rxn 195
  //sp 1
  sp_rates[1] += (fwd_rates[195] - rev_rates[194]);
  //sp 47
  (*dy_N) += (fwd_rates[195] - rev_rates[194]);
  //sp 30
  sp_rates[30] -= (fwd_rates[195] - rev_rates[194]);
  //sp 31
  sp_rates[31] -= (fwd_rates[195] - rev_rates[194]);

  //rxn 196
  //sp 0
  sp_rates[0] += (fwd_rates[196] - rev_rates[195]);
  //sp 5
  sp_rates[5] -= (fwd_rates[196] - rev_rates[195]);
  //sp 38
  sp_rates[38] += (fwd_rates[196] - rev_rates[195]);
  //sp 31
  sp_rates[31] -= (fwd_rates[196] - rev_rates[195]);

  //rxn 197
  //sp 4
  sp_rates[4] += (fwd_rates[197] - rev_rates[196]);
  //sp 35
  sp_rates[35] -= (fwd_rates[197] - rev_rates[196]);
  //sp 47
  (*dy_N) += (fwd_rates[197] - rev_rates[196]);
  //sp 31
  sp_rates[31] -= (fwd_rates[197] - rev_rates[196]);

  //rxn 198
  //sp 1
  sp_rates[1] += (fwd_rates[198] - rev_rates[197]);
  //sp 35
  sp_rates[35] -= (fwd_rates[198] - rev_rates[197]);
  //sp 37
  sp_rates[37] += (fwd_rates[198] - rev_rates[197]);
  //sp 31
  sp_rates[31] -= (fwd_rates[198] - rev_rates[197]);

  //rxn 199
  //sp 32
  sp_rates[32] = -(fwd_rates[199] - rev_rates[198]);
  //sp 2
  sp_rates[2] -= (fwd_rates[199] - rev_rates[198]);
  //sp 4
  sp_rates[4] += (fwd_rates[199] - rev_rates[198]);
  //sp 31
  sp_rates[31] += (fwd_rates[199] - rev_rates[198]);

  //rxn 200
  //sp 32
  sp_rates[32] -= (fwd_rates[200] - rev_rates[199]);
  //sp 1
  sp_rates[1] += (fwd_rates[200] - rev_rates[199]);
  //sp 2
  sp_rates[2] -= (fwd_rates[200] - rev_rates[199]);
  //sp 38
  sp_rates[38] += (fwd_rates[200] - rev_rates[199]);

  //rxn 201
  //sp 32
  sp_rates[32] -= (fwd_rates[201] - rev_rates[200]);
  //sp 1
  sp_rates[1] -= (fwd_rates[201] - rev_rates[200]);
  //sp 0
  sp_rates[0] += (fwd_rates[201] - rev_rates[200]);
  //sp 31
  sp_rates[31] += (fwd_rates[201] - rev_rates[200]);

  //rxn 202
  //sp 32
  sp_rates[32] -= (fwd_rates[202] - rev_rates[201]);
  //sp 4
  sp_rates[4] -= (fwd_rates[202] - rev_rates[201]);
  //sp 5
  sp_rates[5] += (fwd_rates[202] - rev_rates[201]);
  //sp 31
  sp_rates[31] += (fwd_rates[202] - rev_rates[201]);

  //rxn 203
  //sp 1
  sp_rates[1] += (fwd_rates[203] - rev_rates[202]);
  //sp 34
  sp_rates[34] = -(fwd_rates[203] - rev_rates[202]);
  //sp 47
  (*dy_N) += (fwd_rates[203] - rev_rates[202]);

  //rxn 204
  //sp 1
  sp_rates[1] += (fwd_rates[204] - rev_rates[203]) * pres_mod[29];
  //sp 34
  sp_rates[34] -= (fwd_rates[204] - rev_rates[203]) * pres_mod[29];
  //sp 47
  (*dy_N) += (fwd_rates[204] - rev_rates[203]) * pres_mod[29];

  //rxn 205
  //sp 34
  sp_rates[34] -= (fwd_rates[205] - rev_rates[204]);
  //sp 3
  sp_rates[3] -= (fwd_rates[205] - rev_rates[204]);
  //sp 47
  (*dy_N) += (fwd_rates[205] - rev_rates[204]);
  //sp 6
  sp_rates[6] += (fwd_rates[205] - rev_rates[204]);

  //rxn 206
  //sp 4
  sp_rates[4] += (fwd_rates[206] - rev_rates[205]);
  //sp 2
  sp_rates[2] -= (fwd_rates[206] - rev_rates[205]);
  //sp 47
  (*dy_N) += (fwd_rates[206] - rev_rates[205]);
  //sp 34
  sp_rates[34] -= (fwd_rates[206] - rev_rates[205]);

  //rxn 207
  //sp 2
  sp_rates[2] -= (fwd_rates[207] - rev_rates[206]);
  //sp 35
  sp_rates[35] += (fwd_rates[207] - rev_rates[206]);
  //sp 34
  sp_rates[34] -= (fwd_rates[207] - rev_rates[206]);
  //sp 31
  sp_rates[31] += (fwd_rates[207] - rev_rates[206]);

  //rxn 208
  //sp 0
  sp_rates[0] += (fwd_rates[208] - rev_rates[207]);
  //sp 1
  sp_rates[1] -= (fwd_rates[208] - rev_rates[207]);
  //sp 34
  sp_rates[34] -= (fwd_rates[208] - rev_rates[207]);
  //sp 47
  (*dy_N) += (fwd_rates[208] - rev_rates[207]);

  //rxn 209
  //sp 47
  (*dy_N) += (fwd_rates[209] - rev_rates[208]);
  //sp 34
  sp_rates[34] -= (fwd_rates[209] - rev_rates[208]);
  //sp 4
  sp_rates[4] -= (fwd_rates[209] - rev_rates[208]);
  //sp 5
  sp_rates[5] += (fwd_rates[209] - rev_rates[208]);

  //rxn 210
  //sp 47
  (*dy_N) += (fwd_rates[210] - rev_rates[209]);
  //sp 34
  sp_rates[34] -= (fwd_rates[210] - rev_rates[209]);
  //sp 12
  sp_rates[12] -= (fwd_rates[210] - rev_rates[209]);
  //sp 13
  sp_rates[13] += (fwd_rates[210] - rev_rates[209]);

  //rxn 211
  //sp 1
  sp_rates[1] -= (fwd_rates[211] - rev_rates[210]) * pres_mod[30];
  //sp 35
  sp_rates[35] -= (fwd_rates[211] - rev_rates[210]) * pres_mod[30];
  //sp 38
  sp_rates[38] += (fwd_rates[211] - rev_rates[210]) * pres_mod[30];

  //rxn 212
  //sp 2
  sp_rates[2] -= (fwd_rates[212] - rev_rates[211]);
  //sp 35
  sp_rates[35] += (fwd_rates[212] - rev_rates[211]);
  //sp 4
  sp_rates[4] += (fwd_rates[212] - rev_rates[211]);
  //sp 38
  sp_rates[38] -= (fwd_rates[212] - rev_rates[211]);

  //rxn 213
  //sp 0
  sp_rates[0] += (fwd_rates[213] - rev_rates[212]);
  //sp 1
  sp_rates[1] -= (fwd_rates[213] - rev_rates[212]);
  //sp 35
  sp_rates[35] += (fwd_rates[213] - rev_rates[212]);
  //sp 38
  sp_rates[38] -= (fwd_rates[213] - rev_rates[212]);

  //rxn 214
  //sp 35
  sp_rates[35] += (fwd_rates[214] - rev_rates[213]);
  //sp 4
  sp_rates[4] -= (fwd_rates[214] - rev_rates[213]);
  //sp 5
  sp_rates[5] += (fwd_rates[214] - rev_rates[213]);
  //sp 38
  sp_rates[38] -= (fwd_rates[214] - rev_rates[213]);

  //rxn 215
  //sp 35
  sp_rates[35] += (fwd_rates[215] - rev_rates[214]);
  //sp 3
  sp_rates[3] -= (fwd_rates[215] - rev_rates[214]);
  //sp 6
  sp_rates[6] += (fwd_rates[215] - rev_rates[214]);
  //sp 38
  sp_rates[38] -= (fwd_rates[215] - rev_rates[214]);

  //rxn 216
  //sp 2
  sp_rates[2] -= (fwd_rates[216] - rev_rates[215]);
  //sp 30
  sp_rates[30] += (fwd_rates[216] - rev_rates[215]);
  //sp 14
  sp_rates[14] += (fwd_rates[216] - rev_rates[215]);
  //sp 39
  sp_rates[39] = -(fwd_rates[216] - rev_rates[215]);

  //rxn 217
  //sp 1
  sp_rates[1] += (fwd_rates[217] - rev_rates[216]);
  //sp 4
  sp_rates[4] -= (fwd_rates[217] - rev_rates[216]);
  //sp 46
  sp_rates[46] = (fwd_rates[217] - rev_rates[216]);
  //sp 39
  sp_rates[39] -= (fwd_rates[217] - rev_rates[216]);

  //rxn 218
  //sp 40
  sp_rates[40] = (fwd_rates[218] - rev_rates[217]);
  //sp 4
  sp_rates[4] += (fwd_rates[218] - rev_rates[217]);
  //sp 5
  sp_rates[5] -= (fwd_rates[218] - rev_rates[217]);
  //sp 39
  sp_rates[39] -= (fwd_rates[218] - rev_rates[217]);

  //rxn 219
  //sp 2
  sp_rates[2] += (fwd_rates[219] - rev_rates[218]);
  //sp 3
  sp_rates[3] -= (fwd_rates[219] - rev_rates[218]);
  //sp 46
  sp_rates[46] += (fwd_rates[219] - rev_rates[218]);
  //sp 39
  sp_rates[39] -= (fwd_rates[219] - rev_rates[218]);

  //rxn 220
  //sp 0
  sp_rates[0] -= (fwd_rates[220] - rev_rates[219]);
  //sp 1
  sp_rates[1] += (fwd_rates[220] - rev_rates[219]);
  //sp 40
  sp_rates[40] += (fwd_rates[220] - rev_rates[219]);
  //sp 39
  sp_rates[39] -= (fwd_rates[220] - rev_rates[219]);

  //rxn 221
  //sp 2
  sp_rates[2] -= (fwd_rates[221] - rev_rates[220]);
  //sp 35
  sp_rates[35] += (fwd_rates[221] - rev_rates[220]);
  //sp 14
  sp_rates[14] += (fwd_rates[221] - rev_rates[220]);
  //sp 46
  sp_rates[46] -= (fwd_rates[221] - rev_rates[220]);

  //rxn 222
  //sp 1
  sp_rates[1] -= (fwd_rates[222] - rev_rates[221]);
  //sp 14
  sp_rates[14] += (fwd_rates[222] - rev_rates[221]);
  //sp 46
  sp_rates[46] -= (fwd_rates[222] - rev_rates[221]);
  //sp 31
  sp_rates[31] += (fwd_rates[222] - rev_rates[221]);

  //rxn 223
  //sp 1
  sp_rates[1] += (fwd_rates[223] - rev_rates[222]);
  //sp 35
  sp_rates[35] += (fwd_rates[223] - rev_rates[222]);
  //sp 4
  sp_rates[4] -= (fwd_rates[223] - rev_rates[222]);
  //sp 14
  sp_rates[14] += (fwd_rates[223] - rev_rates[222]);
  //sp 46
  sp_rates[46] -= (fwd_rates[223] - rev_rates[222]);

  //rxn 224
  //sp 14
  sp_rates[14] += (fwd_rates[224] - rev_rates[223]);
  //sp 47
  (*dy_N) += (fwd_rates[224] - rev_rates[223]);
  //sp 46
  sp_rates[46] -= (fwd_rates[224] - rev_rates[223]);
  //sp 30
  sp_rates[30] -= (fwd_rates[224] - rev_rates[223]);

  //rxn 225
  //sp 35
  sp_rates[35] += (fwd_rates[225] - rev_rates[224]);
  //sp 3
  sp_rates[3] -= (fwd_rates[225] - rev_rates[224]);
  //sp 46
  sp_rates[46] -= (fwd_rates[225] - rev_rates[224]);
  //sp 15
  sp_rates[15] += (fwd_rates[225] - rev_rates[224]);

  //rxn 226
  //sp 30
  sp_rates[30] += (fwd_rates[226] - rev_rates[225]) * pres_mod[31];
  //sp 14
  sp_rates[14] += (fwd_rates[226] - rev_rates[225]) * pres_mod[31];
  //sp 46
  sp_rates[46] -= (fwd_rates[226] - rev_rates[225]) * pres_mod[31];

  //rxn 227
  //sp 35
  sp_rates[35] -= (fwd_rates[227] - rev_rates[226]);
  //sp 14
  sp_rates[14] += (fwd_rates[227] - rev_rates[226]);
  //sp 46
  sp_rates[46] -= (fwd_rates[227] - rev_rates[226]);
  //sp 37
  sp_rates[37] += (fwd_rates[227] - rev_rates[226]);

  //rxn 228
  //sp 35
  sp_rates[35] -= (fwd_rates[228] - rev_rates[227]);
  //sp 47
  (*dy_N) += (fwd_rates[228] - rev_rates[227]);
  //sp 46
  sp_rates[46] -= (fwd_rates[228] - rev_rates[227]);
  //sp 15
  sp_rates[15] += (fwd_rates[228] - rev_rates[227]);

  //rxn 229
  //sp 40
  sp_rates[40] -= (fwd_rates[229] - rev_rates[228]) * pres_mod[32];
  //sp 1
  sp_rates[1] += (fwd_rates[229] - rev_rates[228]) * pres_mod[32];
  //sp 39
  sp_rates[39] += (fwd_rates[229] - rev_rates[228]) * pres_mod[32];

  //rxn 230
  //sp 40
  sp_rates[40] -= (fwd_rates[230] - rev_rates[229]);
  //sp 1
  sp_rates[1] += (fwd_rates[230] - rev_rates[229]);
  //sp 2
  sp_rates[2] -= (fwd_rates[230] - rev_rates[229]);
  //sp 46
  sp_rates[46] += (fwd_rates[230] - rev_rates[229]);

  //rxn 231
  //sp 40
  sp_rates[40] -= (fwd_rates[231] - rev_rates[230]);
  //sp 2
  sp_rates[2] -= (fwd_rates[231] - rev_rates[230]);
  //sp 14
  sp_rates[14] += (fwd_rates[231] - rev_rates[230]);
  //sp 31
  sp_rates[31] += (fwd_rates[231] - rev_rates[230]);

  //rxn 232
  //sp 40
  sp_rates[40] -= (fwd_rates[232] - rev_rates[231]);
  //sp 2
  sp_rates[2] -= (fwd_rates[232] - rev_rates[231]);
  //sp 4
  sp_rates[4] += (fwd_rates[232] - rev_rates[231]);
  //sp 39
  sp_rates[39] += (fwd_rates[232] - rev_rates[231]);

  //rxn 233
  //sp 40
  sp_rates[40] -= (fwd_rates[233] - rev_rates[232]);
  //sp 44
  sp_rates[44] = (fwd_rates[233] - rev_rates[232]);
  //sp 4
  sp_rates[4] -= (fwd_rates[233] - rev_rates[232]);
  //sp 1
  sp_rates[1] += (fwd_rates[233] - rev_rates[232]);

  //rxn 234
  //sp 40
  sp_rates[40] -= (fwd_rates[234] - rev_rates[233]);
  //sp 1
  sp_rates[1] += (fwd_rates[234] - rev_rates[233]);
  //sp 4
  sp_rates[4] -= (fwd_rates[234] - rev_rates[233]);
  //sp 45
  sp_rates[45] = (fwd_rates[234] - rev_rates[233]);

  //rxn 235
  //sp 40
  sp_rates[40] -= (fwd_rates[235] - rev_rates[234]);
  //sp 32
  sp_rates[32] += (fwd_rates[235] - rev_rates[234]);
  //sp 4
  sp_rates[4] -= (fwd_rates[235] - rev_rates[234]);
  //sp 14
  sp_rates[14] += (fwd_rates[235] - rev_rates[234]);

  //rxn 236
  //sp 40
  sp_rates[40] -= (fwd_rates[236] - rev_rates[235]) * pres_mod[33];
  //sp 1
  sp_rates[1] -= (fwd_rates[236] - rev_rates[235]) * pres_mod[33];
  //sp 41
  sp_rates[41] = (fwd_rates[236] - rev_rates[235]) * pres_mod[33];

  //rxn 237
  //sp 41
  sp_rates[41] -= (fwd_rates[237] - rev_rates[236]);
  //sp 10
  sp_rates[10] += (fwd_rates[237] - rev_rates[236]);
  //sp 47
  (*dy_N) += (fwd_rates[237] - rev_rates[236]);
  //sp 30
  sp_rates[30] -= (fwd_rates[237] - rev_rates[236]);

  //rxn 238
  //sp 8
  sp_rates[8] -= (fwd_rates[238] - rev_rates[237]);
  //sp 47
  (*dy_N) -= (fwd_rates[238] - rev_rates[237]);
  //sp 30
  sp_rates[30] += (fwd_rates[238] - rev_rates[237]);
  //sp 39
  sp_rates[39] += (fwd_rates[238] - rev_rates[237]);

  //rxn 239
  //sp 40
  sp_rates[40] += (fwd_rates[239] - rev_rates[238]);
  //sp 9
  sp_rates[9] -= (fwd_rates[239] - rev_rates[238]);
  //sp 47
  (*dy_N) -= (fwd_rates[239] - rev_rates[238]);
  //sp 30
  sp_rates[30] += (fwd_rates[239] - rev_rates[238]);

  //rxn 240
  //sp 9
  sp_rates[9] -= (fwd_rates[240] - rev_rates[239]) * pres_mod[34];
  //sp 42
  sp_rates[42] = (fwd_rates[240] - rev_rates[239]) * pres_mod[34];
  //sp 47
  (*dy_N) -= (fwd_rates[240] - rev_rates[239]) * pres_mod[34];

  //rxn 241
  //sp 40
  sp_rates[40] += (fwd_rates[241] - rev_rates[240]);
  //sp 10
  sp_rates[10] -= (fwd_rates[241] - rev_rates[240]);
  //sp 47
  (*dy_N) -= (fwd_rates[241] - rev_rates[240]);
  //sp 31
  sp_rates[31] += (fwd_rates[241] - rev_rates[240]);

  //rxn 242
  //sp 40
  sp_rates[40] += (fwd_rates[242] - rev_rates[241]);
  //sp 11
  sp_rates[11] -= (fwd_rates[242] - rev_rates[241]);
  //sp 47
  (*dy_N) -= (fwd_rates[242] - rev_rates[241]);
  //sp 31
  sp_rates[31] += (fwd_rates[242] - rev_rates[241]);

  //rxn 243
  //sp 8
  sp_rates[8] -= (fwd_rates[243] - rev_rates[242]);
  //sp 2
  sp_rates[2] += (fwd_rates[243] - rev_rates[242]);
  //sp 35
  sp_rates[35] -= (fwd_rates[243] - rev_rates[242]);
  //sp 39
  sp_rates[39] += (fwd_rates[243] - rev_rates[242]);

  //rxn 244
  //sp 8
  sp_rates[8] -= (fwd_rates[244] - rev_rates[243]);
  //sp 35
  sp_rates[35] -= (fwd_rates[244] - rev_rates[243]);
  //sp 30
  sp_rates[30] += (fwd_rates[244] - rev_rates[243]);
  //sp 14
  sp_rates[14] += (fwd_rates[244] - rev_rates[243]);

  //rxn 245
  //sp 40
  sp_rates[40] += (fwd_rates[245] - rev_rates[244]);
  //sp 9
  sp_rates[9] -= (fwd_rates[245] - rev_rates[244]);
  //sp 2
  sp_rates[2] += (fwd_rates[245] - rev_rates[244]);
  //sp 35
  sp_rates[35] -= (fwd_rates[245] - rev_rates[244]);

  //rxn 246
  //sp 9
  sp_rates[9] -= (fwd_rates[246] - rev_rates[245]);
  //sp 35
  sp_rates[35] -= (fwd_rates[246] - rev_rates[245]);
  //sp 46
  sp_rates[46] += (fwd_rates[246] - rev_rates[245]);
  //sp 1
  sp_rates[1] += (fwd_rates[246] - rev_rates[245]);

  //rxn 247
  //sp 16
  sp_rates[16] += (fwd_rates[247] - rev_rates[246]);
  //sp 9
  sp_rates[9] -= (fwd_rates[247] - rev_rates[246]);
  //sp 35
  sp_rates[35] -= (fwd_rates[247] - rev_rates[246]);
  //sp 30
  sp_rates[30] += (fwd_rates[247] - rev_rates[246]);

  //rxn 248
  //sp 1
  sp_rates[1] += (fwd_rates[248] - rev_rates[247]);
  //sp 10
  sp_rates[10] -= (fwd_rates[248] - rev_rates[247]);
  //sp 35
  sp_rates[35] -= (fwd_rates[248] - rev_rates[247]);
  //sp 45
  sp_rates[45] += (fwd_rates[248] - rev_rates[247]);

  //rxn 249
  //sp 40
  sp_rates[40] += (fwd_rates[249] - rev_rates[248]);
  //sp 10
  sp_rates[10] -= (fwd_rates[249] - rev_rates[248]);
  //sp 35
  sp_rates[35] -= (fwd_rates[249] - rev_rates[248]);
  //sp 4
  sp_rates[4] += (fwd_rates[249] - rev_rates[248]);

  //rxn 250
  //sp 1
  sp_rates[1] += (fwd_rates[250] - rev_rates[249]);
  //sp 10
  sp_rates[10] -= (fwd_rates[250] - rev_rates[249]);
  //sp 35
  sp_rates[35] -= (fwd_rates[250] - rev_rates[249]);
  //sp 43
  sp_rates[43] = (fwd_rates[250] - rev_rates[249]);

  //rxn 251
  //sp 35
  sp_rates[35] -= (fwd_rates[251] - rev_rates[250]);
  //sp 11
  sp_rates[11] -= (fwd_rates[251] - rev_rates[250]);
  //sp 45
  sp_rates[45] += (fwd_rates[251] - rev_rates[250]);
  //sp 1
  sp_rates[1] += (fwd_rates[251] - rev_rates[250]);

  //rxn 252
  //sp 40
  sp_rates[40] += (fwd_rates[252] - rev_rates[251]);
  //sp 35
  sp_rates[35] -= (fwd_rates[252] - rev_rates[251]);
  //sp 11
  sp_rates[11] -= (fwd_rates[252] - rev_rates[251]);
  //sp 4
  sp_rates[4] += (fwd_rates[252] - rev_rates[251]);

  //rxn 253
  //sp 35
  sp_rates[35] -= (fwd_rates[253] - rev_rates[252]);
  //sp 11
  sp_rates[11] -= (fwd_rates[253] - rev_rates[252]);
  //sp 43
  sp_rates[43] += (fwd_rates[253] - rev_rates[252]);
  //sp 1
  sp_rates[1] += (fwd_rates[253] - rev_rates[252]);

  //rxn 254
  //sp 40
  sp_rates[40] += (fwd_rates[254] - rev_rates[253]);
  //sp 35
  sp_rates[35] -= (fwd_rates[254] - rev_rates[253]);
  //sp 12
  sp_rates[12] -= (fwd_rates[254] - rev_rates[253]);
  //sp 5
  sp_rates[5] += (fwd_rates[254] - rev_rates[253]);

  //rxn 255
  //sp 41
  sp_rates[41] += (fwd_rates[255] - rev_rates[254]);
  //sp 35
  sp_rates[35] -= (fwd_rates[255] - rev_rates[254]);
  //sp 12
  sp_rates[12] -= (fwd_rates[255] - rev_rates[254]);
  //sp 4
  sp_rates[4] += (fwd_rates[255] - rev_rates[254]);

  //rxn 256
  //sp 1
  sp_rates[1] += (fwd_rates[256] - rev_rates[255]);
  //sp 2
  sp_rates[2] -= (fwd_rates[256] - rev_rates[255]);
  //sp 47
  (*dy_N) += (fwd_rates[256] - rev_rates[255]);
  //sp 42
  sp_rates[42] -= (fwd_rates[256] - rev_rates[255]);
  //sp 14
  sp_rates[14] += (fwd_rates[256] - rev_rates[255]);

  //rxn 257
  //sp 40
  sp_rates[40] += (fwd_rates[257] - rev_rates[256]);
  //sp 2
  sp_rates[2] -= (fwd_rates[257] - rev_rates[256]);
  //sp 35
  sp_rates[35] += (fwd_rates[257] - rev_rates[256]);
  //sp 42
  sp_rates[42] -= (fwd_rates[257] - rev_rates[256]);

  //rxn 258
  //sp 16
  sp_rates[16] += (fwd_rates[258] - rev_rates[257]);
  //sp 42
  sp_rates[42] -= (fwd_rates[258] - rev_rates[257]);
  //sp 3
  sp_rates[3] -= (fwd_rates[258] - rev_rates[257]);
  //sp 47
  (*dy_N) += (fwd_rates[258] - rev_rates[257]);
  //sp 2
  sp_rates[2] += (fwd_rates[258] - rev_rates[257]);

  //rxn 259
  //sp 16
  sp_rates[16] += (fwd_rates[259] - rev_rates[258]);
  //sp 1
  sp_rates[1] += (fwd_rates[259] - rev_rates[258]);
  //sp 42
  sp_rates[42] -= (fwd_rates[259] - rev_rates[258]);
  //sp 4
  sp_rates[4] -= (fwd_rates[259] - rev_rates[258]);
  //sp 47
  (*dy_N) += (fwd_rates[259] - rev_rates[258]);

  //rxn 260
  //sp 1
  sp_rates[1] -= (fwd_rates[260] - rev_rates[259]);
  //sp 42
  sp_rates[42] -= (fwd_rates[260] - rev_rates[259]);
  //sp 47
  (*dy_N) += (fwd_rates[260] - rev_rates[259]);
  //sp 10
  sp_rates[10] += (fwd_rates[260] - rev_rates[259]);

  //rxn 261
  //sp 2
  sp_rates[2] -= (fwd_rates[261] - rev_rates[260]);
  //sp 15
  sp_rates[15] += (fwd_rates[261] - rev_rates[260]);
  //sp 45
  sp_rates[45] -= (fwd_rates[261] - rev_rates[260]);
  //sp 31
  sp_rates[31] += (fwd_rates[261] - rev_rates[260]);

  //rxn 262
  //sp 2
  sp_rates[2] -= (fwd_rates[262] - rev_rates[261]);
  //sp 38
  sp_rates[38] += (fwd_rates[262] - rev_rates[261]);
  //sp 45
  sp_rates[45] -= (fwd_rates[262] - rev_rates[261]);
  //sp 14
  sp_rates[14] += (fwd_rates[262] - rev_rates[261]);

  //rxn 263
  //sp 2
  sp_rates[2] -= (fwd_rates[263] - rev_rates[262]);
  //sp 4
  sp_rates[4] += (fwd_rates[263] - rev_rates[262]);
  //sp 45
  sp_rates[45] -= (fwd_rates[263] - rev_rates[262]);
  //sp 46
  sp_rates[46] += (fwd_rates[263] - rev_rates[262]);

  //rxn 264
  //sp 32
  sp_rates[32] += (fwd_rates[264] - rev_rates[263]);
  //sp 1
  sp_rates[1] -= (fwd_rates[264] - rev_rates[263]);
  //sp 45
  sp_rates[45] -= (fwd_rates[264] - rev_rates[263]);
  //sp 14
  sp_rates[14] += (fwd_rates[264] - rev_rates[263]);

  //rxn 265
  //sp 0
  sp_rates[0] += (fwd_rates[265] - rev_rates[264]);
  //sp 1
  sp_rates[1] -= (fwd_rates[265] - rev_rates[264]);
  //sp 45
  sp_rates[45] -= (fwd_rates[265] - rev_rates[264]);
  //sp 46
  sp_rates[46] += (fwd_rates[265] - rev_rates[264]);

  //rxn 266
  //sp 4
  sp_rates[4] -= (fwd_rates[266] - rev_rates[265]);
  //sp 45
  sp_rates[45] -= (fwd_rates[266] - rev_rates[265]);
  //sp 46
  sp_rates[46] += (fwd_rates[266] - rev_rates[265]);
  //sp 5
  sp_rates[5] += (fwd_rates[266] - rev_rates[265]);

  //rxn 267
  //sp 32
  sp_rates[32] += (fwd_rates[267] - rev_rates[266]);
  //sp 4
  sp_rates[4] -= (fwd_rates[267] - rev_rates[266]);
  //sp 45
  sp_rates[45] -= (fwd_rates[267] - rev_rates[266]);
  //sp 15
  sp_rates[15] += (fwd_rates[267] - rev_rates[266]);

  //rxn 268
  //sp 45
  sp_rates[45] -= (fwd_rates[268] - rev_rates[267]) * pres_mod[35];
  //sp 14
  sp_rates[14] += (fwd_rates[268] - rev_rates[267]) * pres_mod[35];
  //sp 31
  sp_rates[31] += (fwd_rates[268] - rev_rates[267]) * pres_mod[35];

  //rxn 269
  //sp 43
  sp_rates[43] -= (fwd_rates[269] - rev_rates[268]);
  //sp 45
  sp_rates[45] += (fwd_rates[269] - rev_rates[268]);

  //rxn 270
  //sp 40
  sp_rates[40] += (fwd_rates[270] - rev_rates[269]);
  //sp 1
  sp_rates[1] -= (fwd_rates[270] - rev_rates[269]);
  //sp 43
  sp_rates[43] -= (fwd_rates[270] - rev_rates[269]);
  //sp 4
  sp_rates[4] += (fwd_rates[270] - rev_rates[269]);

  //rxn 271
  //sp 32
  sp_rates[32] += (fwd_rates[271] - rev_rates[270]);
  //sp 1
  sp_rates[1] -= (fwd_rates[271] - rev_rates[270]);
  //sp 43
  sp_rates[43] -= (fwd_rates[271] - rev_rates[270]);
  //sp 14
  sp_rates[14] += (fwd_rates[271] - rev_rates[270]);

  //rxn 272
  //sp 44
  sp_rates[44] -= (fwd_rates[272] - rev_rates[271]);
  //sp 45
  sp_rates[45] += (fwd_rates[272] - rev_rates[271]);

  //rxn 273
  //sp 35
  sp_rates[35] -= (fwd_rates[273] - rev_rates[272]);
  //sp 27
  sp_rates[27] -= (fwd_rates[273] - rev_rates[272]);
  //sp 14
  sp_rates[14] += (fwd_rates[273] - rev_rates[272]);
  //sp 43
  sp_rates[43] += (fwd_rates[273] - rev_rates[272]);

  //rxn 274
  //sp 1
  sp_rates[1] += (fwd_rates[274] - rev_rates[273]);
  //sp 12
  sp_rates[12] -= (fwd_rates[274] - rev_rates[273]);
  //sp 30
  sp_rates[30] -= (fwd_rates[274] - rev_rates[273]);
  //sp 41
  sp_rates[41] += (fwd_rates[274] - rev_rates[273]);

  //rxn 275
  //sp 0
  sp_rates[0] += (fwd_rates[275] - rev_rates[274]);
  //sp 40
  sp_rates[40] += (fwd_rates[275] - rev_rates[274]);
  //sp 12
  sp_rates[12] -= (fwd_rates[275] - rev_rates[274]);
  //sp 30
  sp_rates[30] -= (fwd_rates[275] - rev_rates[274]);

  //rxn 276
  //sp 0
  sp_rates[0] += (fwd_rates[276] - rev_rates[275]);
  //sp 1
  sp_rates[1] -= (fwd_rates[276] - rev_rates[275]);
  //sp 32
  sp_rates[32] += (fwd_rates[276] - rev_rates[275]);
  //sp 33
  sp_rates[33] = -(fwd_rates[276] - rev_rates[275]);

  //rxn 277
  //sp 32
  sp_rates[32] += (fwd_rates[277] - rev_rates[276]);
  //sp 33
  sp_rates[33] -= (fwd_rates[277] - rev_rates[276]);
  //sp 4
  sp_rates[4] -= (fwd_rates[277] - rev_rates[276]);
  //sp 5
  sp_rates[5] += (fwd_rates[277] - rev_rates[276]);

  //rxn 278
  //sp 32
  sp_rates[32] += (fwd_rates[278] - rev_rates[277]);
  //sp 33
  sp_rates[33] -= (fwd_rates[278] - rev_rates[277]);
  //sp 2
  sp_rates[2] -= (fwd_rates[278] - rev_rates[277]);
  //sp 4
  sp_rates[4] += (fwd_rates[278] - rev_rates[277]);

  //rxn 279
  //sp 15
  sp_rates[15] -= (fwd_rates[279] - rev_rates[278]);
  //sp 38
  sp_rates[38] += (fwd_rates[279] - rev_rates[278]);
  //sp 14
  sp_rates[14] += (fwd_rates[279] - rev_rates[278]);
  //sp 31
  sp_rates[31] -= (fwd_rates[279] - rev_rates[278]);

  //rxn 280
  //sp 35
  sp_rates[35] += (fwd_rates[280] - rev_rates[279]);
  //sp 36
  sp_rates[36] -= (fwd_rates[280] - rev_rates[279]);
  //sp 46
  sp_rates[46] += (fwd_rates[280] - rev_rates[279]);
  //sp 39
  sp_rates[39] -= (fwd_rates[280] - rev_rates[279]);

  //rxn 281
  //sp 36
  sp_rates[36] -= (fwd_rates[281] - rev_rates[280]);
  //sp 37
  sp_rates[37] += (fwd_rates[281] - rev_rates[280]);
  //sp 46
  sp_rates[46] -= (fwd_rates[281] - rev_rates[280]);
  //sp 15
  sp_rates[15] += (fwd_rates[281] - rev_rates[280]);

  //rxn 282
  //sp 35
  sp_rates[35] += (fwd_rates[282] - rev_rates[281]);
  //sp 14
  sp_rates[14] += (fwd_rates[282] - rev_rates[281]);
  //sp 30
  sp_rates[30] -= (fwd_rates[282] - rev_rates[281]);
  //sp 15
  sp_rates[15] -= (fwd_rates[282] - rev_rates[281]);

  //rxn 283
  //sp 0
  sp_rates[0] += fwd_rates[283];
  //sp 1
  sp_rates[1] += fwd_rates[283];
  //sp 2
  sp_rates[2] -= fwd_rates[283];
  //sp 12
  sp_rates[12] -= fwd_rates[283];
  //sp 14
  sp_rates[14] += fwd_rates[283];

  //rxn 284
  //sp 24
  sp_rates[24] -= (fwd_rates[284] - rev_rates[282]);
  //sp 1
  sp_rates[1] += (fwd_rates[284] - rev_rates[282]);
  //sp 2
  sp_rates[2] -= (fwd_rates[284] - rev_rates[282]);
  //sp 51
  sp_rates[50] = (fwd_rates[284] - rev_rates[282]);

  //rxn 285
  //sp 25
  sp_rates[25] -= (fwd_rates[285] - rev_rates[283]);
  //sp 2
  sp_rates[2] -= (fwd_rates[285] - rev_rates[283]);
  //sp 52
  sp_rates[51] = (fwd_rates[285] - rev_rates[283]);
  //sp 1
  sp_rates[1] += (fwd_rates[285] - rev_rates[283]);

  //rxn 286
  //sp 3
  sp_rates[3] += (fwd_rates[286] - rev_rates[284]);
  //sp 4
  sp_rates[4] -= (fwd_rates[286] - rev_rates[284]);
  //sp 5
  sp_rates[5] += (fwd_rates[286] - rev_rates[284]);
  //sp 6
  sp_rates[6] -= (fwd_rates[286] - rev_rates[284]);

  //rxn 287
  //sp 0
  sp_rates[0] += fwd_rates[287];
  //sp 4
  sp_rates[4] -= fwd_rates[287];
  //sp 12
  sp_rates[12] -= fwd_rates[287];
  //sp 17
  sp_rates[17] += fwd_rates[287];

  //rxn 288
  //sp 0
  sp_rates[0] -= (fwd_rates[288] - rev_rates[285]) * pres_mod[36];
  //sp 9
  sp_rates[9] -= (fwd_rates[288] - rev_rates[285]) * pres_mod[36];
  //sp 12
  sp_rates[12] += (fwd_rates[288] - rev_rates[285]) * pres_mod[36];

  //rxn 289
  //sp 1
  sp_rates[1] += 2.0 * fwd_rates[289];
  //sp 10
  sp_rates[10] -= fwd_rates[289];
  //sp 3
  sp_rates[3] -= fwd_rates[289];
  //sp 15
  sp_rates[15] += fwd_rates[289];

  //rxn 290
  //sp 17
  sp_rates[17] += (fwd_rates[290] - rev_rates[286]);
  //sp 10
  sp_rates[10] -= (fwd_rates[290] - rev_rates[286]);
  //sp 3
  sp_rates[3] -= (fwd_rates[290] - rev_rates[286]);
  //sp 2
  sp_rates[2] += (fwd_rates[290] - rev_rates[286]);

  //rxn 291
  //sp 1
  sp_rates[1] += 2.0 * fwd_rates[291];
  //sp 10
  sp_rates[10] -= 2.0 * fwd_rates[291];
  //sp 22
  sp_rates[22] += fwd_rates[291];

  //rxn 292
  //sp 0
  sp_rates[0] += fwd_rates[292];
  //sp 17
  sp_rates[17] += fwd_rates[292];
  //sp 11
  sp_rates[11] -= fwd_rates[292];
  //sp 5
  sp_rates[5] -= fwd_rates[292];

  //rxn 293
  //sp 51
  sp_rates[50] += (fwd_rates[293] - rev_rates[287]);
  //sp 3
  sp_rates[3] -= (fwd_rates[293] - rev_rates[287]);
  //sp 2
  sp_rates[2] += (fwd_rates[293] - rev_rates[287]);
  //sp 23
  sp_rates[23] -= (fwd_rates[293] - rev_rates[287]);

  //rxn 294
  //sp 3
  sp_rates[3] -= (fwd_rates[294] - rev_rates[288]);
  //sp 22
  sp_rates[22] += (fwd_rates[294] - rev_rates[288]);
  //sp 6
  sp_rates[6] += (fwd_rates[294] - rev_rates[288]);
  //sp 23
  sp_rates[23] -= (fwd_rates[294] - rev_rates[288]);

  //rxn 295
  //sp 2
  sp_rates[2] -= (fwd_rates[295] - rev_rates[289]);
  //sp 52
  sp_rates[51] -= (fwd_rates[295] - rev_rates[289]);
  //sp 4
  sp_rates[4] += (fwd_rates[295] - rev_rates[289]);
  //sp 51
  sp_rates[50] += (fwd_rates[295] - rev_rates[289]);

  //rxn 296
  //sp 4
  sp_rates[4] += fwd_rates[296];
  //sp 2
  sp_rates[2] -= fwd_rates[296];
  //sp 52
  sp_rates[51] -= fwd_rates[296];
  //sp 12
  sp_rates[12] += fwd_rates[296];
  //sp 14
  sp_rates[14] += fwd_rates[296];

  //rxn 297
  //sp 3
  sp_rates[3] -= fwd_rates[297];
  //sp 52
  sp_rates[51] -= fwd_rates[297];
  //sp 12
  sp_rates[12] += fwd_rates[297];
  //sp 6
  sp_rates[6] += fwd_rates[297];
  //sp 14
  sp_rates[14] += fwd_rates[297];

  //rxn 298
  //sp 0
  sp_rates[0] += (fwd_rates[298] - rev_rates[290]);
  //sp 1
  sp_rates[1] -= (fwd_rates[298] - rev_rates[290]);
  //sp 51
  sp_rates[50] += (fwd_rates[298] - rev_rates[290]);
  //sp 52
  sp_rates[51] -= (fwd_rates[298] - rev_rates[290]);

  //rxn 299
  //sp 0
  sp_rates[0] += fwd_rates[299];
  //sp 1
  sp_rates[1] -= fwd_rates[299];
  //sp 52
  sp_rates[51] -= fwd_rates[299];
  //sp 12
  sp_rates[12] += fwd_rates[299];
  //sp 14
  sp_rates[14] += fwd_rates[299];

  //rxn 300
  //sp 12
  sp_rates[12] += fwd_rates[300];
  //sp 52
  sp_rates[51] -= fwd_rates[300];
  //sp 4
  sp_rates[4] -= fwd_rates[300];
  //sp 5
  sp_rates[5] += fwd_rates[300];
  //sp 14
  sp_rates[14] += fwd_rates[300];

  //rxn 301
  //sp 52
  sp_rates[51] -= fwd_rates[301];
  //sp 12
  sp_rates[12] += fwd_rates[301];
  //sp 14
  sp_rates[14] += fwd_rates[301];
  //sp 6
  sp_rates[6] -= fwd_rates[301];
  //sp 7
  sp_rates[7] += fwd_rates[301];

  //rxn 302
  //sp 52
  sp_rates[51] -= fwd_rates[302];
  //sp 13
  sp_rates[13] += fwd_rates[302];
  //sp 14
  sp_rates[14] += fwd_rates[302];

  //rxn 303
  //sp 1
  sp_rates[1] -= (fwd_rates[303] - rev_rates[291]) * pres_mod[37];
  //sp 51
  sp_rates[50] += (fwd_rates[303] - rev_rates[291]) * pres_mod[37];
  //sp 28
  sp_rates[28] -= (fwd_rates[303] - rev_rates[291]) * pres_mod[37];

  //rxn 304
  //sp 1
  sp_rates[1] += fwd_rates[304];
  //sp 51
  sp_rates[50] -= fwd_rates[304];
  //sp 10
  sp_rates[10] += fwd_rates[304];
  //sp 2
  sp_rates[2] -= fwd_rates[304];
  //sp 15
  sp_rates[15] += fwd_rates[304];

  //rxn 305
  //sp 17
  sp_rates[17] += fwd_rates[305];
  //sp 51
  sp_rates[50] -= fwd_rates[305];
  //sp 3
  sp_rates[3] -= fwd_rates[305];
  //sp 4
  sp_rates[4] += fwd_rates[305];
  //sp 14
  sp_rates[14] += fwd_rates[305];

  //rxn 306
  //sp 16
  sp_rates[16] += 2.0 * fwd_rates[306];
  //sp 51
  sp_rates[50] -= fwd_rates[306];
  //sp 3
  sp_rates[3] -= fwd_rates[306];
  //sp 4
  sp_rates[4] += fwd_rates[306];

  //rxn 307
  //sp 16
  sp_rates[16] += (fwd_rates[307] - rev_rates[292]);
  //sp 1
  sp_rates[1] -= (fwd_rates[307] - rev_rates[292]);
  //sp 51
  sp_rates[50] -= (fwd_rates[307] - rev_rates[292]);
  //sp 12
  sp_rates[12] += (fwd_rates[307] - rev_rates[292]);

  //rxn 308
  //sp 0
  sp_rates[0] += (fwd_rates[308] - rev_rates[293]);
  //sp 1
  sp_rates[1] -= (fwd_rates[308] - rev_rates[293]);
  //sp 51
  sp_rates[50] -= (fwd_rates[308] - rev_rates[293]);
  //sp 28
  sp_rates[28] += (fwd_rates[308] - rev_rates[293]);

  //rxn 309
  //sp 28
  sp_rates[28] += (fwd_rates[309] - rev_rates[294]);
  //sp 51
  sp_rates[50] -= (fwd_rates[309] - rev_rates[294]);
  //sp 4
  sp_rates[4] -= (fwd_rates[309] - rev_rates[294]);
  //sp 5
  sp_rates[5] += (fwd_rates[309] - rev_rates[294]);

  //rxn 310
  //sp 16
  sp_rates[16] += (fwd_rates[310] - rev_rates[295]);
  //sp 51
  sp_rates[50] -= (fwd_rates[310] - rev_rates[295]);
  //sp 4
  sp_rates[4] -= (fwd_rates[310] - rev_rates[295]);
  //sp 18
  sp_rates[18] += (fwd_rates[310] - rev_rates[295]);

  //rxn 311
  //sp 25
  sp_rates[25] -= (fwd_rates[311] - rev_rates[296]) * pres_mod[38];
  //sp 12
  sp_rates[12] -= (fwd_rates[311] - rev_rates[296]) * pres_mod[38];
  //sp 50
  sp_rates[49] = (fwd_rates[311] - rev_rates[296]) * pres_mod[38];

  //rxn 312
  //sp 49
  sp_rates[48] = (fwd_rates[312] - rev_rates[297]);
  //sp 50
  sp_rates[49] -= (fwd_rates[312] - rev_rates[297]);
  //sp 2
  sp_rates[2] -= (fwd_rates[312] - rev_rates[297]);
  //sp 4
  sp_rates[4] += (fwd_rates[312] - rev_rates[297]);

  //rxn 313
  //sp 0
  sp_rates[0] += (fwd_rates[313] - rev_rates[298]);
  //sp 1
  sp_rates[1] -= (fwd_rates[313] - rev_rates[298]);
  //sp 49
  sp_rates[48] += (fwd_rates[313] - rev_rates[298]);
  //sp 50
  sp_rates[49] -= (fwd_rates[313] - rev_rates[298]);

  //rxn 314
  //sp 49
  sp_rates[48] += (fwd_rates[314] - rev_rates[299]);
  //sp 50
  sp_rates[49] -= (fwd_rates[314] - rev_rates[299]);
  //sp 4
  sp_rates[4] -= (fwd_rates[314] - rev_rates[299]);
  //sp 5
  sp_rates[5] += (fwd_rates[314] - rev_rates[299]);

  //rxn 315
  //sp 49
  sp_rates[48] -= (fwd_rates[315] - rev_rates[300]);
  //sp 50
  sp_rates[49] += (fwd_rates[315] - rev_rates[300]);
  //sp 6
  sp_rates[6] += (fwd_rates[315] - rev_rates[300]);
  //sp 7
  sp_rates[7] -= (fwd_rates[315] - rev_rates[300]);

  //rxn 316
  //sp 49
  sp_rates[48] += (fwd_rates[316] - rev_rates[301]);
  //sp 50
  sp_rates[49] -= (fwd_rates[316] - rev_rates[301]);
  //sp 12
  sp_rates[12] -= (fwd_rates[316] - rev_rates[301]);
  //sp 13
  sp_rates[13] += (fwd_rates[316] - rev_rates[301]);

  //rxn 317
  //sp 24
  sp_rates[24] -= (fwd_rates[317] - rev_rates[302]) * pres_mod[39];
  //sp 49
  sp_rates[48] += (fwd_rates[317] - rev_rates[302]) * pres_mod[39];
  //sp 12
  sp_rates[12] -= (fwd_rates[317] - rev_rates[302]) * pres_mod[39];

  //rxn 318
  //sp 49
  sp_rates[48] -= (fwd_rates[318] - rev_rates[303]);
  //sp 17
  sp_rates[17] += (fwd_rates[318] - rev_rates[303]);
  //sp 2
  sp_rates[2] -= (fwd_rates[318] - rev_rates[303]);
  //sp 25
  sp_rates[25] += (fwd_rates[318] - rev_rates[303]);

  //rxn 319
  //sp 49
  sp_rates[48] -= (fwd_rates[319] - rev_rates[304]) * pres_mod[40];
  //sp 1
  sp_rates[1] -= (fwd_rates[319] - rev_rates[304]) * pres_mod[40];
  //sp 50
  sp_rates[49] += (fwd_rates[319] - rev_rates[304]) * pres_mod[40];

  //rxn 320
  //sp 49
  sp_rates[48] -= (fwd_rates[320] - rev_rates[305]);
  //sp 1
  sp_rates[1] -= (fwd_rates[320] - rev_rates[305]);
  //sp 12
  sp_rates[12] += (fwd_rates[320] - rev_rates[305]);
  //sp 25
  sp_rates[25] += (fwd_rates[320] - rev_rates[305]);

  //rxn 321
  //sp 49
  sp_rates[48] -= (fwd_rates[321] - rev_rates[306]);
  //sp 25
  sp_rates[25] += (fwd_rates[321] - rev_rates[306]);
  //sp 18
  sp_rates[18] += (fwd_rates[321] - rev_rates[306]);
  //sp 4
  sp_rates[4] -= (fwd_rates[321] - rev_rates[306]);

  //rxn 322
  //sp 49
  sp_rates[48] -= (fwd_rates[322] - rev_rates[307]);
  //sp 50
  sp_rates[49] += (fwd_rates[322] - rev_rates[307]);
  //sp 3
  sp_rates[3] += (fwd_rates[322] - rev_rates[307]);
  //sp 6
  sp_rates[6] -= (fwd_rates[322] - rev_rates[307]);

  //rxn 323
  //sp 49
  sp_rates[48] -= fwd_rates[323];
  //sp 17
  sp_rates[17] += fwd_rates[323];
  //sp 4
  sp_rates[4] += fwd_rates[323];
  //sp 6
  sp_rates[6] -= fwd_rates[323];
  //sp 25
  sp_rates[25] += fwd_rates[323];

  //rxn 324
  //sp 49
  sp_rates[48] -= (fwd_rates[324] - rev_rates[308]);
  //sp 25
  sp_rates[25] += 2.0 * (fwd_rates[324] - rev_rates[308]);
  //sp 12
  sp_rates[12] -= (fwd_rates[324] - rev_rates[308]);

  //sp 48
  sp_rates[47] = 0.0;
} // end eval_spec_rates

