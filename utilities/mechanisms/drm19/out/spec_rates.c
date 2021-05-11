#include "header.h"
#include "rates.h"

void eval_spec_rates (const double * __restrict__ fwd_rates, const double * __restrict__ rev_rates, const double * __restrict__ pres_mod, double * __restrict__ sp_rates, double * __restrict__ dy_N) {
  //rxn 0
  //sp 1
  sp_rates[1] = -(fwd_rates[0] - rev_rates[0]) * pres_mod[0];
  //sp 2
  sp_rates[2] = -(fwd_rates[0] - rev_rates[0]) * pres_mod[0];
  //sp 4
  sp_rates[4] = (fwd_rates[0] - rev_rates[0]) * pres_mod[0];

  //rxn 1
  //sp 0
  sp_rates[0] = -(fwd_rates[1] - rev_rates[1]);
  //sp 1
  sp_rates[1] += (fwd_rates[1] - rev_rates[1]);
  //sp 2
  sp_rates[2] -= (fwd_rates[1] - rev_rates[1]);
  //sp 4
  sp_rates[4] += (fwd_rates[1] - rev_rates[1]);

  //rxn 2
  //sp 2
  sp_rates[2] -= (fwd_rates[2] - rev_rates[2]);
  //sp 3
  sp_rates[3] = (fwd_rates[2] - rev_rates[2]);
  //sp 4
  sp_rates[4] += (fwd_rates[2] - rev_rates[2]);
  //sp 6
  sp_rates[6] = -(fwd_rates[2] - rev_rates[2]);

  //rxn 3
  //sp 1
  sp_rates[1] += (fwd_rates[3] - rev_rates[3]);
  //sp 2
  sp_rates[2] -= (fwd_rates[3] - rev_rates[3]);
  //sp 13
  sp_rates[13] = (fwd_rates[3] - rev_rates[3]);
  //sp 7
  sp_rates[7] = -(fwd_rates[3] - rev_rates[3]);

  //rxn 4
  //sp 8
  sp_rates[8] = -(fwd_rates[4] - rev_rates[4]);
  //sp 1
  sp_rates[1] += (fwd_rates[4] - rev_rates[4]);
  //sp 2
  sp_rates[2] -= (fwd_rates[4] - rev_rates[4]);
  //sp 13
  sp_rates[13] += (fwd_rates[4] - rev_rates[4]);

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
  //sp 9
  sp_rates[9] += (fwd_rates[6] - rev_rates[6]);
  //sp 10
  sp_rates[10] = -(fwd_rates[6] - rev_rates[6]);
  //sp 4
  sp_rates[4] += (fwd_rates[6] - rev_rates[6]);
  //sp 2
  sp_rates[2] -= (fwd_rates[6] - rev_rates[6]);

  //rxn 7
  //sp 2
  sp_rates[2] -= (fwd_rates[7] - rev_rates[7]) * pres_mod[1];
  //sp 11
  sp_rates[11] = -(fwd_rates[7] - rev_rates[7]) * pres_mod[1];
  //sp 12
  sp_rates[12] = (fwd_rates[7] - rev_rates[7]) * pres_mod[1];

  //rxn 8
  //sp 2
  sp_rates[2] -= (fwd_rates[8] - rev_rates[8]);
  //sp 11
  sp_rates[11] += (fwd_rates[8] - rev_rates[8]);
  //sp 4
  sp_rates[4] += (fwd_rates[8] - rev_rates[8]);
  //sp 13
  sp_rates[13] -= (fwd_rates[8] - rev_rates[8]);

  //rxn 9
  //sp 1
  sp_rates[1] += (fwd_rates[9] - rev_rates[9]);
  //sp 2
  sp_rates[2] -= (fwd_rates[9] - rev_rates[9]);
  //sp 12
  sp_rates[12] += (fwd_rates[9] - rev_rates[9]);
  //sp 13
  sp_rates[13] -= (fwd_rates[9] - rev_rates[9]);

  //rxn 10
  //sp 2
  sp_rates[2] -= (fwd_rates[10] - rev_rates[10]);
  //sp 4
  sp_rates[4] += (fwd_rates[10] - rev_rates[10]);
  //sp 13
  sp_rates[13] += (fwd_rates[10] - rev_rates[10]);
  //sp 14
  sp_rates[14] -= (fwd_rates[10] - rev_rates[10]);

  //rxn 11
  //sp 16
  sp_rates[16] = -(fwd_rates[11] - rev_rates[11]);
  //sp 9
  sp_rates[9] += (fwd_rates[11] - rev_rates[11]);
  //sp 2
  sp_rates[2] -= (fwd_rates[11] - rev_rates[11]);
  //sp 13
  sp_rates[13] += (fwd_rates[11] - rev_rates[11]);

  //rxn 12
  //sp 17
  sp_rates[17] = -(fwd_rates[12] - rev_rates[12]);
  //sp 2
  sp_rates[2] -= (fwd_rates[12] - rev_rates[12]);
  //sp 14
  sp_rates[14] += (fwd_rates[12] - rev_rates[12]);
  //sp 9
  sp_rates[9] += (fwd_rates[12] - rev_rates[12]);

  //rxn 13
  //sp 17
  sp_rates[17] += (fwd_rates[13] - rev_rates[13]);
  //sp 18
  sp_rates[18] = -(fwd_rates[13] - rev_rates[13]);
  //sp 4
  sp_rates[4] += (fwd_rates[13] - rev_rates[13]);
  //sp 2
  sp_rates[2] -= (fwd_rates[13] - rev_rates[13]);

  //rxn 14
  //sp 3
  sp_rates[3] -= (fwd_rates[14] - rev_rates[14]);
  //sp 2
  sp_rates[2] += (fwd_rates[14] - rev_rates[14]);
  //sp 11
  sp_rates[11] -= (fwd_rates[14] - rev_rates[14]);
  //sp 12
  sp_rates[12] += (fwd_rates[14] - rev_rates[14]);

  //rxn 15
  //sp 3
  sp_rates[3] -= (fwd_rates[15] - rev_rates[15]);
  //sp 6
  sp_rates[6] += (fwd_rates[15] - rev_rates[15]);
  //sp 14
  sp_rates[14] -= (fwd_rates[15] - rev_rates[15]);
  //sp 13
  sp_rates[13] += (fwd_rates[15] - rev_rates[15]);

  //rxn 16
  //sp 1
  sp_rates[1] -= (fwd_rates[16] - rev_rates[16]) * pres_mod[2];
  //sp 3
  sp_rates[3] -= (fwd_rates[16] - rev_rates[16]) * pres_mod[2];
  //sp 6
  sp_rates[6] += (fwd_rates[16] - rev_rates[16]) * pres_mod[2];

  //rxn 17
  //sp 1
  sp_rates[1] -= (fwd_rates[17] - rev_rates[17]);
  //sp 3
  sp_rates[3] -= (fwd_rates[17] - rev_rates[17]);
  //sp 6
  sp_rates[6] += (fwd_rates[17] - rev_rates[17]);

  //rxn 18
  //sp 1
  sp_rates[1] -= (fwd_rates[18] - rev_rates[18]);
  //sp 3
  sp_rates[3] -= (fwd_rates[18] - rev_rates[18]);
  //sp 6
  sp_rates[6] += (fwd_rates[18] - rev_rates[18]);

  //rxn 19
  //sp 1
  sp_rates[1] -= (fwd_rates[19] - rev_rates[19]);
  //sp 3
  sp_rates[3] -= (fwd_rates[19] - rev_rates[19]);
  //sp 6
  sp_rates[6] += (fwd_rates[19] - rev_rates[19]);

  //rxn 20
  //sp 3
  sp_rates[3] -= (fwd_rates[20] - rev_rates[20]);
  //sp 1
  sp_rates[1] -= (fwd_rates[20] - rev_rates[20]);
  //sp 6
  sp_rates[6] += (fwd_rates[20] - rev_rates[20]);

  //rxn 21
  //sp 1
  sp_rates[1] -= (fwd_rates[21] - rev_rates[21]);
  //sp 2
  sp_rates[2] += (fwd_rates[21] - rev_rates[21]);
  //sp 3
  sp_rates[3] -= (fwd_rates[21] - rev_rates[21]);
  //sp 4
  sp_rates[4] += (fwd_rates[21] - rev_rates[21]);

  //rxn 22
  //sp 0
  sp_rates[0] += (fwd_rates[22] - rev_rates[22]) * pres_mod[3];
  //sp 1
  sp_rates[1] -= 2.0 * (fwd_rates[22] - rev_rates[22]) * pres_mod[3];

  //rxn 23
  //sp 0
  sp_rates[0] += (fwd_rates[23] - rev_rates[23]);
  //sp 1
  sp_rates[1] -= 2.0 * (fwd_rates[23] - rev_rates[23]);

  //rxn 24
  //sp 0
  sp_rates[0] += (fwd_rates[24] - rev_rates[24]);
  //sp 1
  sp_rates[1] -= 2.0 * (fwd_rates[24] - rev_rates[24]);

  //rxn 25
  //sp 0
  sp_rates[0] += (fwd_rates[25] - rev_rates[25]);
  //sp 1
  sp_rates[1] -= 2.0 * (fwd_rates[25] - rev_rates[25]);

  //rxn 26
  //sp 1
  sp_rates[1] -= (fwd_rates[26] - rev_rates[26]) * pres_mod[4];
  //sp 4
  sp_rates[4] -= (fwd_rates[26] - rev_rates[26]) * pres_mod[4];
  //sp 5
  sp_rates[5] = (fwd_rates[26] - rev_rates[26]) * pres_mod[4];

  //rxn 27
  //sp 0
  sp_rates[0] += (fwd_rates[27] - rev_rates[27]);
  //sp 1
  sp_rates[1] -= (fwd_rates[27] - rev_rates[27]);
  //sp 3
  sp_rates[3] += (fwd_rates[27] - rev_rates[27]);
  //sp 6
  sp_rates[6] -= (fwd_rates[27] - rev_rates[27]);

  //rxn 28
  //sp 1
  sp_rates[1] -= (fwd_rates[28] - rev_rates[28]);
  //sp 4
  sp_rates[4] += 2.0 * (fwd_rates[28] - rev_rates[28]);
  //sp 6
  sp_rates[6] -= (fwd_rates[28] - rev_rates[28]);

  //rxn 29
  //sp 1
  sp_rates[1] -= (fwd_rates[29] - rev_rates[29]) * pres_mod[5];
  //sp 9
  sp_rates[9] += (fwd_rates[29] - rev_rates[29]) * pres_mod[5];
  //sp 7
  sp_rates[7] -= (fwd_rates[29] - rev_rates[29]) * pres_mod[5];

  //rxn 30
  //sp 1
  sp_rates[1] -= (fwd_rates[30] - rev_rates[30]) * pres_mod[6];
  //sp 10
  sp_rates[10] += (fwd_rates[30] - rev_rates[30]) * pres_mod[6];
  //sp 9
  sp_rates[9] -= (fwd_rates[30] - rev_rates[30]) * pres_mod[6];

  //rxn 31
  //sp 0
  sp_rates[0] += (fwd_rates[31] - rev_rates[31]);
  //sp 1
  sp_rates[1] -= (fwd_rates[31] - rev_rates[31]);
  //sp 10
  sp_rates[10] -= (fwd_rates[31] - rev_rates[31]);
  //sp 9
  sp_rates[9] += (fwd_rates[31] - rev_rates[31]);

  //rxn 32
  //sp 1
  sp_rates[1] -= (fwd_rates[32] - rev_rates[32]) * pres_mod[7];
  //sp 13
  sp_rates[13] -= (fwd_rates[32] - rev_rates[32]) * pres_mod[7];
  //sp 14
  sp_rates[14] += (fwd_rates[32] - rev_rates[32]) * pres_mod[7];

  //rxn 33
  //sp 0
  sp_rates[0] += (fwd_rates[33] - rev_rates[33]);
  //sp 1
  sp_rates[1] -= (fwd_rates[33] - rev_rates[33]);
  //sp 11
  sp_rates[11] += (fwd_rates[33] - rev_rates[33]);
  //sp 13
  sp_rates[13] -= (fwd_rates[33] - rev_rates[33]);

  //rxn 34
  //sp 1
  sp_rates[1] -= (fwd_rates[34] - rev_rates[34]) * pres_mod[8];
  //sp 14
  sp_rates[14] -= (fwd_rates[34] - rev_rates[34]) * pres_mod[8];
  //sp 15
  sp_rates[15] = (fwd_rates[34] - rev_rates[34]) * pres_mod[8];

  //rxn 35
  //sp 0
  sp_rates[0] += (fwd_rates[35] - rev_rates[35]);
  //sp 1
  sp_rates[1] -= (fwd_rates[35] - rev_rates[35]);
  //sp 13
  sp_rates[13] += (fwd_rates[35] - rev_rates[35]);
  //sp 14
  sp_rates[14] -= (fwd_rates[35] - rev_rates[35]);

  //rxn 36
  //sp 1
  sp_rates[1] -= (fwd_rates[36] - rev_rates[36]);
  //sp 9
  sp_rates[9] += (fwd_rates[36] - rev_rates[36]);
  //sp 4
  sp_rates[4] += (fwd_rates[36] - rev_rates[36]);
  //sp 15
  sp_rates[15] -= (fwd_rates[36] - rev_rates[36]);

  //rxn 37
  //sp 16
  sp_rates[16] -= (fwd_rates[37] - rev_rates[37]) * pres_mod[9];
  //sp 1
  sp_rates[1] -= (fwd_rates[37] - rev_rates[37]) * pres_mod[9];
  //sp 17
  sp_rates[17] += (fwd_rates[37] - rev_rates[37]) * pres_mod[9];

  //rxn 38
  //sp 1
  sp_rates[1] -= (fwd_rates[38] - rev_rates[38]) * pres_mod[10];
  //sp 18
  sp_rates[18] += (fwd_rates[38] - rev_rates[38]) * pres_mod[10];
  //sp 17
  sp_rates[17] -= (fwd_rates[38] - rev_rates[38]) * pres_mod[10];

  //rxn 39
  //sp 0
  sp_rates[0] += (fwd_rates[39] - rev_rates[39]);
  //sp 1
  sp_rates[1] -= (fwd_rates[39] - rev_rates[39]);
  //sp 18
  sp_rates[18] -= (fwd_rates[39] - rev_rates[39]);
  //sp 17
  sp_rates[17] += (fwd_rates[39] - rev_rates[39]);

  //rxn 40
  //sp 0
  sp_rates[0] -= (fwd_rates[40] - rev_rates[40]) * pres_mod[11];
  //sp 11
  sp_rates[11] -= (fwd_rates[40] - rev_rates[40]) * pres_mod[11];
  //sp 14
  sp_rates[14] += (fwd_rates[40] - rev_rates[40]) * pres_mod[11];

  //rxn 41
  //sp 0
  sp_rates[0] -= (fwd_rates[41] - rev_rates[41]);
  //sp 1
  sp_rates[1] += (fwd_rates[41] - rev_rates[41]);
  //sp 4
  sp_rates[4] -= (fwd_rates[41] - rev_rates[41]);
  //sp 5
  sp_rates[5] += (fwd_rates[41] - rev_rates[41]);

  //rxn 42
  //sp 2
  sp_rates[2] += (fwd_rates[42] - rev_rates[42]);
  //sp 4
  sp_rates[4] -= 2.0 * (fwd_rates[42] - rev_rates[42]);
  //sp 5
  sp_rates[5] += (fwd_rates[42] - rev_rates[42]);

  //rxn 43
  //sp 3
  sp_rates[3] += (fwd_rates[43] - rev_rates[43]);
  //sp 4
  sp_rates[4] -= (fwd_rates[43] - rev_rates[43]);
  //sp 5
  sp_rates[5] += (fwd_rates[43] - rev_rates[43]);
  //sp 6
  sp_rates[6] -= (fwd_rates[43] - rev_rates[43]);

  //rxn 44
  //sp 1
  sp_rates[1] += (fwd_rates[44] - rev_rates[44]);
  //sp 4
  sp_rates[4] -= (fwd_rates[44] - rev_rates[44]);
  //sp 14
  sp_rates[14] += (fwd_rates[44] - rev_rates[44]);
  //sp 7
  sp_rates[7] -= (fwd_rates[44] - rev_rates[44]);

  //rxn 45
  //sp 8
  sp_rates[8] -= (fwd_rates[45] - rev_rates[45]);
  //sp 1
  sp_rates[1] += (fwd_rates[45] - rev_rates[45]);
  //sp 4
  sp_rates[4] -= (fwd_rates[45] - rev_rates[45]);
  //sp 14
  sp_rates[14] += (fwd_rates[45] - rev_rates[45]);

  //rxn 46
  //sp 9
  sp_rates[9] -= (fwd_rates[46] - rev_rates[46]);
  //sp 4
  sp_rates[4] -= (fwd_rates[46] - rev_rates[46]);
  //sp 5
  sp_rates[5] += (fwd_rates[46] - rev_rates[46]);
  //sp 7
  sp_rates[7] += (fwd_rates[46] - rev_rates[46]);

  //rxn 47
  //sp 8
  sp_rates[8] += (fwd_rates[47] - rev_rates[47]);
  //sp 9
  sp_rates[9] -= (fwd_rates[47] - rev_rates[47]);
  //sp 4
  sp_rates[4] -= (fwd_rates[47] - rev_rates[47]);
  //sp 5
  sp_rates[5] += (fwd_rates[47] - rev_rates[47]);

  //rxn 48
  //sp 9
  sp_rates[9] += (fwd_rates[48] - rev_rates[48]);
  //sp 10
  sp_rates[10] -= (fwd_rates[48] - rev_rates[48]);
  //sp 4
  sp_rates[4] -= (fwd_rates[48] - rev_rates[48]);
  //sp 5
  sp_rates[5] += (fwd_rates[48] - rev_rates[48]);

  //rxn 49
  //sp 1
  sp_rates[1] += (fwd_rates[49] - rev_rates[49]);
  //sp 11
  sp_rates[11] -= (fwd_rates[49] - rev_rates[49]);
  //sp 4
  sp_rates[4] -= (fwd_rates[49] - rev_rates[49]);
  //sp 12
  sp_rates[12] += (fwd_rates[49] - rev_rates[49]);

  //rxn 50
  //sp 11
  sp_rates[11] += (fwd_rates[50] - rev_rates[50]);
  //sp 4
  sp_rates[4] -= (fwd_rates[50] - rev_rates[50]);
  //sp 13
  sp_rates[13] -= (fwd_rates[50] - rev_rates[50]);
  //sp 5
  sp_rates[5] += (fwd_rates[50] - rev_rates[50]);

  //rxn 51
  //sp 4
  sp_rates[4] -= (fwd_rates[51] - rev_rates[51]);
  //sp 5
  sp_rates[5] += (fwd_rates[51] - rev_rates[51]);
  //sp 14
  sp_rates[14] -= (fwd_rates[51] - rev_rates[51]);
  //sp 13
  sp_rates[13] += (fwd_rates[51] - rev_rates[51]);

  //rxn 52
  //sp 17
  sp_rates[17] += (fwd_rates[52] - rev_rates[52]);
  //sp 18
  sp_rates[18] -= (fwd_rates[52] - rev_rates[52]);
  //sp 4
  sp_rates[4] -= (fwd_rates[52] - rev_rates[52]);
  //sp 5
  sp_rates[5] += (fwd_rates[52] - rev_rates[52]);

  //rxn 53
  //sp 4
  sp_rates[4] += (fwd_rates[53] - rev_rates[53]);
  //sp 14
  sp_rates[14] += (fwd_rates[53] - rev_rates[53]);
  //sp 6
  sp_rates[6] -= (fwd_rates[53] - rev_rates[53]);
  //sp 7
  sp_rates[7] -= (fwd_rates[53] - rev_rates[53]);

  //rxn 54
  //sp 9
  sp_rates[9] -= (fwd_rates[54] - rev_rates[54]);
  //sp 10
  sp_rates[10] += (fwd_rates[54] - rev_rates[54]);
  //sp 3
  sp_rates[3] += (fwd_rates[54] - rev_rates[54]);
  //sp 6
  sp_rates[6] -= (fwd_rates[54] - rev_rates[54]);

  //rxn 55
  //sp 9
  sp_rates[9] -= (fwd_rates[55] - rev_rates[55]);
  //sp 4
  sp_rates[4] += (fwd_rates[55] - rev_rates[55]);
  //sp 6
  sp_rates[6] -= (fwd_rates[55] - rev_rates[55]);
  //sp 15
  sp_rates[15] += (fwd_rates[55] - rev_rates[55]);

  //rxn 56
  //sp 4
  sp_rates[4] += (fwd_rates[56] - rev_rates[56]);
  //sp 11
  sp_rates[11] -= (fwd_rates[56] - rev_rates[56]);
  //sp 12
  sp_rates[12] += (fwd_rates[56] - rev_rates[56]);
  //sp 6
  sp_rates[6] -= (fwd_rates[56] - rev_rates[56]);

  //rxn 57
  //sp 3
  sp_rates[3] -= (fwd_rates[57] - rev_rates[57]);
  //sp 4
  sp_rates[4] += (fwd_rates[57] - rev_rates[57]);
  //sp 13
  sp_rates[13] += (fwd_rates[57] - rev_rates[57]);
  //sp 7
  sp_rates[7] -= (fwd_rates[57] - rev_rates[57]);

  //rxn 58
  //sp 0
  sp_rates[0] -= (fwd_rates[58] - rev_rates[58]);
  //sp 1
  sp_rates[1] += (fwd_rates[58] - rev_rates[58]);
  //sp 9
  sp_rates[9] += (fwd_rates[58] - rev_rates[58]);
  //sp 7
  sp_rates[7] -= (fwd_rates[58] - rev_rates[58]);

  //rxn 59
  //sp 16
  sp_rates[16] += (fwd_rates[59] - rev_rates[59]);
  //sp 9
  sp_rates[9] -= (fwd_rates[59] - rev_rates[59]);
  //sp 1
  sp_rates[1] += (fwd_rates[59] - rev_rates[59]);
  //sp 7
  sp_rates[7] -= (fwd_rates[59] - rev_rates[59]);

  //rxn 60
  //sp 9
  sp_rates[9] += 2.0 * (fwd_rates[60] - rev_rates[60]);
  //sp 10
  sp_rates[10] -= (fwd_rates[60] - rev_rates[60]);
  //sp 7
  sp_rates[7] -= (fwd_rates[60] - rev_rates[60]);

  //rxn 61
  //sp 8
  sp_rates[8] -= (fwd_rates[61] - rev_rates[61]);
  //sp 7
  sp_rates[7] += (fwd_rates[61] - rev_rates[61]);

  //rxn 62
  //sp 8
  sp_rates[8] -= (fwd_rates[62] - rev_rates[62]);
  //sp 7
  sp_rates[7] += (fwd_rates[62] - rev_rates[62]);

  //rxn 63
  //sp 8
  sp_rates[8] -= (fwd_rates[63] - rev_rates[63]);
  //sp 1
  sp_rates[1] += (fwd_rates[63] - rev_rates[63]);
  //sp 3
  sp_rates[3] -= (fwd_rates[63] - rev_rates[63]);
  //sp 4
  sp_rates[4] += (fwd_rates[63] - rev_rates[63]);
  //sp 11
  sp_rates[11] += (fwd_rates[63] - rev_rates[63]);

  //rxn 64
  //sp 8
  sp_rates[8] -= (fwd_rates[64] - rev_rates[64]);
  //sp 11
  sp_rates[11] += (fwd_rates[64] - rev_rates[64]);
  //sp 3
  sp_rates[3] -= (fwd_rates[64] - rev_rates[64]);
  //sp 5
  sp_rates[5] += (fwd_rates[64] - rev_rates[64]);

  //rxn 65
  //sp 0
  sp_rates[0] -= (fwd_rates[65] - rev_rates[65]);
  //sp 8
  sp_rates[8] -= (fwd_rates[65] - rev_rates[65]);
  //sp 9
  sp_rates[9] += (fwd_rates[65] - rev_rates[65]);
  //sp 1
  sp_rates[1] += (fwd_rates[65] - rev_rates[65]);

  //rxn 66
  //sp 8
  sp_rates[8] -= (fwd_rates[66] - rev_rates[66]);
  //sp 7
  sp_rates[7] += (fwd_rates[66] - rev_rates[66]);

  //rxn 67
  //sp 8
  sp_rates[8] -= (fwd_rates[67] - rev_rates[67]);
  //sp 9
  sp_rates[9] -= (fwd_rates[67] - rev_rates[67]);
  //sp 16
  sp_rates[16] += (fwd_rates[67] - rev_rates[67]);
  //sp 1
  sp_rates[1] += (fwd_rates[67] - rev_rates[67]);

  //rxn 68
  //sp 8
  sp_rates[8] -= (fwd_rates[68] - rev_rates[68]);
  //sp 9
  sp_rates[9] += 2.0 * (fwd_rates[68] - rev_rates[68]);
  //sp 10
  sp_rates[10] -= (fwd_rates[68] - rev_rates[68]);

  //rxn 69
  //sp 8
  sp_rates[8] -= (fwd_rates[69] - rev_rates[69]);
  //sp 7
  sp_rates[7] += (fwd_rates[69] - rev_rates[69]);

  //rxn 70
  //sp 8
  sp_rates[8] -= (fwd_rates[70] - rev_rates[70]);
  //sp 7
  sp_rates[7] += (fwd_rates[70] - rev_rates[70]);

  //rxn 71
  //sp 8
  sp_rates[8] -= (fwd_rates[71] - rev_rates[71]);
  //sp 11
  sp_rates[11] += (fwd_rates[71] - rev_rates[71]);
  //sp 12
  sp_rates[12] -= (fwd_rates[71] - rev_rates[71]);
  //sp 14
  sp_rates[14] += (fwd_rates[71] - rev_rates[71]);

  //rxn 72
  //sp 9
  sp_rates[9] -= (fwd_rates[72] - rev_rates[72]);
  //sp 2
  sp_rates[2] += (fwd_rates[72] - rev_rates[72]);
  //sp 3
  sp_rates[3] -= (fwd_rates[72] - rev_rates[72]);
  //sp 15
  sp_rates[15] += (fwd_rates[72] - rev_rates[72]);

  //rxn 73
  //sp 9
  sp_rates[9] -= (fwd_rates[73] - rev_rates[73]);
  //sp 3
  sp_rates[3] -= (fwd_rates[73] - rev_rates[73]);
  //sp 4
  sp_rates[4] += (fwd_rates[73] - rev_rates[73]);
  //sp 14
  sp_rates[14] += (fwd_rates[73] - rev_rates[73]);

  //rxn 74
  //sp 9
  sp_rates[9] -= 2.0 * (fwd_rates[74] - rev_rates[74]) * pres_mod[12];
  //sp 18
  sp_rates[18] += (fwd_rates[74] - rev_rates[74]) * pres_mod[12];

  //rxn 75
  //sp 9
  sp_rates[9] -= 2.0 * (fwd_rates[75] - rev_rates[75]);
  //sp 17
  sp_rates[17] += (fwd_rates[75] - rev_rates[75]);
  //sp 1
  sp_rates[1] += (fwd_rates[75] - rev_rates[75]);

  //rxn 76
  //sp 9
  sp_rates[9] -= (fwd_rates[76] - rev_rates[76]);
  //sp 10
  sp_rates[10] += (fwd_rates[76] - rev_rates[76]);
  //sp 11
  sp_rates[11] += (fwd_rates[76] - rev_rates[76]);
  //sp 13
  sp_rates[13] -= (fwd_rates[76] - rev_rates[76]);

  //rxn 77
  //sp 9
  sp_rates[9] -= (fwd_rates[77] - rev_rates[77]);
  //sp 10
  sp_rates[10] += (fwd_rates[77] - rev_rates[77]);
  //sp 13
  sp_rates[13] += (fwd_rates[77] - rev_rates[77]);
  //sp 14
  sp_rates[14] -= (fwd_rates[77] - rev_rates[77]);

  //rxn 78
  //sp 9
  sp_rates[9] -= (fwd_rates[78] - rev_rates[78]);
  //sp 18
  sp_rates[18] -= (fwd_rates[78] - rev_rates[78]);
  //sp 10
  sp_rates[10] += (fwd_rates[78] - rev_rates[78]);
  //sp 17
  sp_rates[17] += (fwd_rates[78] - rev_rates[78]);

  //rxn 79
  //sp 1
  sp_rates[1] += (fwd_rates[79] - rev_rates[79]);
  //sp 11
  sp_rates[11] += (fwd_rates[79] - rev_rates[79]);
  //sp 13
  sp_rates[13] -= (fwd_rates[79] - rev_rates[79]);

  //rxn 80
  //sp 1
  sp_rates[1] += (fwd_rates[80] - rev_rates[80]) * pres_mod[13];
  //sp 11
  sp_rates[11] += (fwd_rates[80] - rev_rates[80]) * pres_mod[13];
  //sp 13
  sp_rates[13] -= (fwd_rates[80] - rev_rates[80]) * pres_mod[13];

  //rxn 81
  //sp 11
  sp_rates[11] += (fwd_rates[81] - rev_rates[81]);
  //sp 3
  sp_rates[3] -= (fwd_rates[81] - rev_rates[81]);
  //sp 13
  sp_rates[13] -= (fwd_rates[81] - rev_rates[81]);
  //sp 6
  sp_rates[6] += (fwd_rates[81] - rev_rates[81]);

  //rxn 82
  //sp 3
  sp_rates[3] -= (fwd_rates[82] - rev_rates[82]);
  //sp 6
  sp_rates[6] += (fwd_rates[82] - rev_rates[82]);
  //sp 14
  sp_rates[14] += (fwd_rates[82] - rev_rates[82]);
  //sp 15
  sp_rates[15] -= (fwd_rates[82] - rev_rates[82]);

  //rxn 83
  //sp 16
  sp_rates[16] += (fwd_rates[83] - rev_rates[83]);
  //sp 17
  sp_rates[17] -= (fwd_rates[83] - rev_rates[83]);
  //sp 3
  sp_rates[3] -= (fwd_rates[83] - rev_rates[83]);
  //sp 6
  sp_rates[6] += (fwd_rates[83] - rev_rates[83]);

  //sp 20
  sp_rates[19] = 0.0;
  //sp 19
  (*dy_N) = 0.0;
} // end eval_spec_rates

