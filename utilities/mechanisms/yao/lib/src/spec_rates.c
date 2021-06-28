#include "header.h"
#include "rates.h"

void eval_spec_rates (const double * __restrict__ fwd_rates, const double * __restrict__ rev_rates, const double * __restrict__ pres_mod, double * __restrict__ sp_rates, double * __restrict__ dy_N) {
  //rxn 0
  //sp 9
  sp_rates[8] = -(fwd_rates[0] - rev_rates[0]);
  //sp 2
  sp_rates[1] = -(fwd_rates[0] - rev_rates[0]);
  //sp 3
  sp_rates[2] = (fwd_rates[0] - rev_rates[0]);
  //sp 4
  sp_rates[3] = (fwd_rates[0] - rev_rates[0]);

  //rxn 1
  //sp 2
  sp_rates[1] += (fwd_rates[1] - rev_rates[1]);
  //sp 3
  sp_rates[2] -= (fwd_rates[1] - rev_rates[1]);
  //sp 4
  sp_rates[3] += (fwd_rates[1] - rev_rates[1]);
  //sp 6
  sp_rates[5] = -(fwd_rates[1] - rev_rates[1]);

  //rxn 2
  //sp 2
  sp_rates[1] += (fwd_rates[2] - rev_rates[2]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2] - rev_rates[2]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2] - rev_rates[2]);
  //sp 7
  sp_rates[6] = (fwd_rates[2] - rev_rates[2]);

  //rxn 3
  //sp 3
  sp_rates[2] += (fwd_rates[3] - rev_rates[3]);
  //sp 4
  sp_rates[3] -= 2.0 * (fwd_rates[3] - rev_rates[3]);
  //sp 7
  sp_rates[6] += (fwd_rates[3] - rev_rates[3]);

  //rxn 4
  //sp 9
  sp_rates[8] -= (fwd_rates[4] - rev_rates[4]) * pres_mod[0];
  //sp 2
  sp_rates[1] -= (fwd_rates[4] - rev_rates[4]) * pres_mod[0];
  //sp 5
  sp_rates[4] = (fwd_rates[4] - rev_rates[4]) * pres_mod[0];

  //rxn 5
  //sp 2
  sp_rates[1] -= (fwd_rates[5] - rev_rates[5]);
  //sp 4
  sp_rates[3] += 2.0 * (fwd_rates[5] - rev_rates[5]);
  //sp 5
  sp_rates[4] -= (fwd_rates[5] - rev_rates[5]);

  //rxn 6
  //sp 9
  sp_rates[8] -= (fwd_rates[6] - rev_rates[6]);
  //sp 2
  sp_rates[1] += (fwd_rates[6] - rev_rates[6]);
  //sp 5
  sp_rates[4] += (fwd_rates[6] - rev_rates[6]);
  //sp 6
  sp_rates[5] -= (fwd_rates[6] - rev_rates[6]);

  //rxn 7
  //sp 9
  sp_rates[8] += (fwd_rates[7] - rev_rates[7]);
  //sp 4
  sp_rates[3] -= (fwd_rates[7] - rev_rates[7]);
  //sp 5
  sp_rates[4] -= (fwd_rates[7] - rev_rates[7]);
  //sp 7
  sp_rates[6] += (fwd_rates[7] - rev_rates[7]);

  //rxn 8
  //sp 2
  sp_rates[1] -= (fwd_rates[8] - rev_rates[8]);
  //sp 3
  sp_rates[2] += (fwd_rates[8] - rev_rates[8]);
  //sp 5
  sp_rates[4] -= (fwd_rates[8] - rev_rates[8]);
  //sp 7
  sp_rates[6] += (fwd_rates[8] - rev_rates[8]);

  //rxn 9
  //sp 9
  sp_rates[8] += (fwd_rates[9] - rev_rates[9]);
  //sp 3
  sp_rates[2] -= (fwd_rates[9] - rev_rates[9]);
  //sp 4
  sp_rates[3] += (fwd_rates[9] - rev_rates[9]);
  //sp 5
  sp_rates[4] -= (fwd_rates[9] - rev_rates[9]);

  //rxn 10
  //sp 9
  sp_rates[8] += (fwd_rates[10] - rev_rates[10]);
  //sp 5
  sp_rates[4] -= 2.0 * (fwd_rates[10] - rev_rates[10]);
  //sp 8
  sp_rates[7] = (fwd_rates[10] - rev_rates[10]);

  //rxn 11
  //sp 9
  sp_rates[8] += (fwd_rates[11] - rev_rates[11]);
  //sp 5
  sp_rates[4] -= 2.0 * (fwd_rates[11] - rev_rates[11]);
  //sp 8
  sp_rates[7] += (fwd_rates[11] - rev_rates[11]);

  //rxn 12
  //sp 2
  sp_rates[1] -= (fwd_rates[12] - rev_rates[12]);
  //sp 4
  sp_rates[3] += (fwd_rates[12] - rev_rates[12]);
  //sp 7
  sp_rates[6] += (fwd_rates[12] - rev_rates[12]);
  //sp 8
  sp_rates[7] -= (fwd_rates[12] - rev_rates[12]);

  //rxn 13
  //sp 2
  sp_rates[1] -= (fwd_rates[13] - rev_rates[13]);
  //sp 5
  sp_rates[4] += (fwd_rates[13] - rev_rates[13]);
  //sp 6
  sp_rates[5] += (fwd_rates[13] - rev_rates[13]);
  //sp 8
  sp_rates[7] -= (fwd_rates[13] - rev_rates[13]);

  //rxn 14
  //sp 3
  sp_rates[2] -= (fwd_rates[14] - rev_rates[14]);
  //sp 4
  sp_rates[3] += (fwd_rates[14] - rev_rates[14]);
  //sp 5
  sp_rates[4] += (fwd_rates[14] - rev_rates[14]);
  //sp 8
  sp_rates[7] -= (fwd_rates[14] - rev_rates[14]);

  //rxn 15
  //sp 4
  sp_rates[3] -= (fwd_rates[15] - rev_rates[15]);
  //sp 5
  sp_rates[4] += (fwd_rates[15] - rev_rates[15]);
  //sp 7
  sp_rates[6] += (fwd_rates[15] - rev_rates[15]);
  //sp 8
  sp_rates[7] -= (fwd_rates[15] - rev_rates[15]);

  //rxn 16
  //sp 4
  sp_rates[3] -= (fwd_rates[16] - rev_rates[16]);
  //sp 5
  sp_rates[4] += (fwd_rates[16] - rev_rates[16]);
  //sp 7
  sp_rates[6] += (fwd_rates[16] - rev_rates[16]);
  //sp 8
  sp_rates[7] -= (fwd_rates[16] - rev_rates[16]);

  //rxn 17
  //sp 4
  sp_rates[3] -= 2.0 * (fwd_rates[17] - rev_rates[17]) * pres_mod[1];
  //sp 8
  sp_rates[7] += (fwd_rates[17] - rev_rates[17]) * pres_mod[1];

  //rxn 18
  //sp 2
  sp_rates[1] -= 2.0 * (fwd_rates[18] - rev_rates[18]) * pres_mod[2];
  //sp 6
  sp_rates[5] += (fwd_rates[18] - rev_rates[18]) * pres_mod[2];

  //rxn 19
  //sp 2
  sp_rates[1] -= (fwd_rates[19] - rev_rates[19]) * pres_mod[3];
  //sp 4
  sp_rates[3] -= (fwd_rates[19] - rev_rates[19]) * pres_mod[3];
  //sp 7
  sp_rates[6] += (fwd_rates[19] - rev_rates[19]) * pres_mod[3];

  //rxn 20
  //sp 9
  sp_rates[8] += (fwd_rates[20] - rev_rates[20]) * pres_mod[4];
  //sp 3
  sp_rates[2] -= 2.0 * (fwd_rates[20] - rev_rates[20]) * pres_mod[4];

  //rxn 21
  //sp 2
  sp_rates[1] -= 2.0 * (fwd_rates[21] - rev_rates[21]);
  //sp 6
  sp_rates[5] += (fwd_rates[21] - rev_rates[21]);

  //rxn 22
  //sp 2
  sp_rates[1] -= 2.0 * (fwd_rates[22] - rev_rates[22]);
  //sp 6
  sp_rates[5] += (fwd_rates[22] - rev_rates[22]);

  //rxn 23
  //sp 2
  sp_rates[1] -= 2.0 * (fwd_rates[23] - rev_rates[23]);
  //sp 6
  sp_rates[5] += (fwd_rates[23] - rev_rates[23]);

  //rxn 24
  //sp 2
  sp_rates[1] -= (fwd_rates[24] - rev_rates[24]) * pres_mod[5];
  //sp 3
  sp_rates[2] -= (fwd_rates[24] - rev_rates[24]) * pres_mod[5];
  //sp 4
  sp_rates[3] += (fwd_rates[24] - rev_rates[24]) * pres_mod[5];

  //rxn 25
  //sp 17
  sp_rates[16] = -(fwd_rates[25] - rev_rates[25]);
  //sp 2
  sp_rates[1] += (fwd_rates[25] - rev_rates[25]);
  //sp 4
  sp_rates[3] -= (fwd_rates[25] - rev_rates[25]);
  //sp 18
  sp_rates[17] = (fwd_rates[25] - rev_rates[25]);

  //rxn 26
  //sp 17
  sp_rates[16] -= (fwd_rates[26] - rev_rates[26]);
  //sp 2
  sp_rates[1] += (fwd_rates[26] - rev_rates[26]);
  //sp 4
  sp_rates[3] -= (fwd_rates[26] - rev_rates[26]);
  //sp 18
  sp_rates[17] += (fwd_rates[26] - rev_rates[26]);

  //rxn 27
  //sp 17
  sp_rates[16] -= (fwd_rates[27] - rev_rates[27]);
  //sp 18
  sp_rates[17] += (fwd_rates[27] - rev_rates[27]);
  //sp 4
  sp_rates[3] += (fwd_rates[27] - rev_rates[27]);
  //sp 5
  sp_rates[4] -= (fwd_rates[27] - rev_rates[27]);

  //rxn 28
  //sp 17
  sp_rates[16] -= (fwd_rates[28] - rev_rates[28]) * pres_mod[6];
  //sp 18
  sp_rates[17] += (fwd_rates[28] - rev_rates[28]) * pres_mod[6];
  //sp 3
  sp_rates[2] -= (fwd_rates[28] - rev_rates[28]) * pres_mod[6];

  //rxn 29
  //sp 17
  sp_rates[16] -= (fwd_rates[29] - rev_rates[29]);
  //sp 9
  sp_rates[8] -= (fwd_rates[29] - rev_rates[29]);
  //sp 3
  sp_rates[2] += (fwd_rates[29] - rev_rates[29]);
  //sp 18
  sp_rates[17] += (fwd_rates[29] - rev_rates[29]);

  //rxn 30
  //sp 17
  sp_rates[16] += (fwd_rates[30] - rev_rates[30]) * pres_mod[7];
  //sp 2
  sp_rates[1] += (fwd_rates[30] - rev_rates[30]) * pres_mod[7];
  //sp 14
  sp_rates[13] = -(fwd_rates[30] - rev_rates[30]) * pres_mod[7];

  //rxn 31
  //sp 17
  sp_rates[16] += (fwd_rates[31] - rev_rates[31]);
  //sp 2
  sp_rates[1] -= (fwd_rates[31] - rev_rates[31]);
  //sp 14
  sp_rates[13] -= (fwd_rates[31] - rev_rates[31]);
  //sp 6
  sp_rates[5] += (fwd_rates[31] - rev_rates[31]);

  //rxn 32
  //sp 17
  sp_rates[16] += (fwd_rates[32] - rev_rates[32]);
  //sp 3
  sp_rates[2] -= (fwd_rates[32] - rev_rates[32]);
  //sp 4
  sp_rates[3] += (fwd_rates[32] - rev_rates[32]);
  //sp 14
  sp_rates[13] -= (fwd_rates[32] - rev_rates[32]);

  //rxn 33
  //sp 2
  sp_rates[1] += (fwd_rates[33] - rev_rates[33]);
  //sp 3
  sp_rates[2] -= (fwd_rates[33] - rev_rates[33]);
  //sp 14
  sp_rates[13] -= (fwd_rates[33] - rev_rates[33]);
  //sp 18
  sp_rates[17] += (fwd_rates[33] - rev_rates[33]);

  //rxn 34
  //sp 17
  sp_rates[16] += (fwd_rates[34] - rev_rates[34]);
  //sp 4
  sp_rates[3] -= (fwd_rates[34] - rev_rates[34]);
  //sp 14
  sp_rates[13] -= (fwd_rates[34] - rev_rates[34]);
  //sp 7
  sp_rates[6] += (fwd_rates[34] - rev_rates[34]);

  //rxn 35
  //sp 9
  sp_rates[8] -= (fwd_rates[35] - rev_rates[35]);
  //sp 17
  sp_rates[16] += (fwd_rates[35] - rev_rates[35]);
  //sp 5
  sp_rates[4] += (fwd_rates[35] - rev_rates[35]);
  //sp 14
  sp_rates[13] -= (fwd_rates[35] - rev_rates[35]);

  //rxn 36
  //sp 17
  sp_rates[16] += (fwd_rates[36] - rev_rates[36]);
  //sp 2
  sp_rates[1] += (fwd_rates[36] - rev_rates[36]);
  //sp 14
  sp_rates[13] -= (fwd_rates[36] - rev_rates[36]);

  //rxn 37
  //sp 17
  sp_rates[16] -= (fwd_rates[37] - rev_rates[37]) * pres_mod[8];
  //sp 6
  sp_rates[5] -= (fwd_rates[37] - rev_rates[37]) * pres_mod[8];
  //sp 15
  sp_rates[14] = (fwd_rates[37] - rev_rates[37]) * pres_mod[8];

  //rxn 38
  //sp 2
  sp_rates[1] -= (fwd_rates[38] - rev_rates[38]) * pres_mod[9];
  //sp 14
  sp_rates[13] -= (fwd_rates[38] - rev_rates[38]) * pres_mod[9];
  //sp 15
  sp_rates[14] += (fwd_rates[38] - rev_rates[38]) * pres_mod[9];

  //rxn 39
  //sp 2
  sp_rates[1] -= (fwd_rates[39] - rev_rates[39]) * pres_mod[10];
  //sp 12
  sp_rates[11] = (fwd_rates[39] - rev_rates[39]) * pres_mod[10];
  //sp 10
  sp_rates[9] = -(fwd_rates[39] - rev_rates[39]) * pres_mod[10];

  //rxn 40
  //sp 10
  sp_rates[9] -= (fwd_rates[40] - rev_rates[40]);
  //sp 3
  sp_rates[2] -= (fwd_rates[40] - rev_rates[40]);
  //sp 14
  sp_rates[13] += (fwd_rates[40] - rev_rates[40]);
  //sp 2
  sp_rates[1] += (fwd_rates[40] - rev_rates[40]);

  //rxn 41
  //sp 10
  sp_rates[9] -= (fwd_rates[41] - rev_rates[41]);
  //sp 4
  sp_rates[3] -= (fwd_rates[41] - rev_rates[41]);
  //sp 15
  sp_rates[14] += (fwd_rates[41] - rev_rates[41]);
  //sp 2
  sp_rates[1] += (fwd_rates[41] - rev_rates[41]);

  //rxn 42
  //sp 10
  sp_rates[9] -= (fwd_rates[42] - rev_rates[42]);
  //sp 12
  sp_rates[11] += (fwd_rates[42] - rev_rates[42]);
  //sp 6
  sp_rates[5] -= (fwd_rates[42] - rev_rates[42]);
  //sp 2
  sp_rates[1] += (fwd_rates[42] - rev_rates[42]);

  //rxn 43
  //sp 9
  sp_rates[8] -= (fwd_rates[43] - rev_rates[43]);
  //sp 10
  sp_rates[9] -= (fwd_rates[43] - rev_rates[43]);
  //sp 4
  sp_rates[3] += (fwd_rates[43] - rev_rates[43]);
  //sp 14
  sp_rates[13] += (fwd_rates[43] - rev_rates[43]);

  //rxn 44
  //sp 9
  sp_rates[8] -= (fwd_rates[44] - rev_rates[44]);
  //sp 10
  sp_rates[9] -= (fwd_rates[44] - rev_rates[44]);
  //sp 18
  sp_rates[17] += (fwd_rates[44] - rev_rates[44]);
  //sp 2
  sp_rates[1] += 2.0 * (fwd_rates[44] - rev_rates[44]);

  //rxn 45
  //sp 10
  sp_rates[9] -= (fwd_rates[45] - rev_rates[45]);
  //sp 4
  sp_rates[3] += (fwd_rates[45] - rev_rates[45]);
  //sp 5
  sp_rates[4] -= (fwd_rates[45] - rev_rates[45]);
  //sp 15
  sp_rates[14] += (fwd_rates[45] - rev_rates[45]);

  //rxn 46
  //sp 10
  sp_rates[9] -= 2.0 * (fwd_rates[46] - rev_rates[46]);
  //sp 19
  sp_rates[18] = (fwd_rates[46] - rev_rates[46]);
  //sp 6
  sp_rates[5] += (fwd_rates[46] - rev_rates[46]);

  //rxn 47
  //sp 10
  sp_rates[9] += (fwd_rates[47] - rev_rates[47]);
  //sp 11
  sp_rates[10] = -(fwd_rates[47] - rev_rates[47]);

  //rxn 48
  //sp 10
  sp_rates[9] += (fwd_rates[48] - rev_rates[48]);
  //sp 11
  sp_rates[10] -= (fwd_rates[48] - rev_rates[48]);

  //rxn 49
  //sp 17
  sp_rates[16] += (fwd_rates[49] - rev_rates[49]);
  //sp 11
  sp_rates[10] -= (fwd_rates[49] - rev_rates[49]);
  //sp 3
  sp_rates[2] -= (fwd_rates[49] - rev_rates[49]);
  //sp 6
  sp_rates[5] += (fwd_rates[49] - rev_rates[49]);

  //rxn 50
  //sp 2
  sp_rates[1] += (fwd_rates[50] - rev_rates[50]);
  //sp 11
  sp_rates[10] -= (fwd_rates[50] - rev_rates[50]);
  //sp 3
  sp_rates[2] -= (fwd_rates[50] - rev_rates[50]);
  //sp 14
  sp_rates[13] += (fwd_rates[50] - rev_rates[50]);

  //rxn 51
  //sp 2
  sp_rates[1] += (fwd_rates[51] - rev_rates[51]);
  //sp 11
  sp_rates[10] -= (fwd_rates[51] - rev_rates[51]);
  //sp 4
  sp_rates[3] -= (fwd_rates[51] - rev_rates[51]);
  //sp 15
  sp_rates[14] += (fwd_rates[51] - rev_rates[51]);

  //rxn 52
  //sp 2
  sp_rates[1] += (fwd_rates[52] - rev_rates[52]);
  //sp 11
  sp_rates[10] -= (fwd_rates[52] - rev_rates[52]);
  //sp 12
  sp_rates[11] += (fwd_rates[52] - rev_rates[52]);
  //sp 6
  sp_rates[5] -= (fwd_rates[52] - rev_rates[52]);

  //rxn 53
  //sp 9
  sp_rates[8] -= (fwd_rates[53] - rev_rates[53]);
  //sp 2
  sp_rates[1] += (fwd_rates[53] - rev_rates[53]);
  //sp 11
  sp_rates[10] -= (fwd_rates[53] - rev_rates[53]);
  //sp 4
  sp_rates[3] += (fwd_rates[53] - rev_rates[53]);
  //sp 17
  sp_rates[16] += (fwd_rates[53] - rev_rates[53]);

  //rxn 54
  //sp 9
  sp_rates[8] -= (fwd_rates[54] - rev_rates[54]);
  //sp 17
  sp_rates[16] += (fwd_rates[54] - rev_rates[54]);
  //sp 11
  sp_rates[10] -= (fwd_rates[54] - rev_rates[54]);
  //sp 7
  sp_rates[6] += (fwd_rates[54] - rev_rates[54]);

  //rxn 55
  //sp 10
  sp_rates[9] += (fwd_rates[55] - rev_rates[55]);
  //sp 11
  sp_rates[10] -= (fwd_rates[55] - rev_rates[55]);

  //rxn 56
  //sp 10
  sp_rates[9] += (fwd_rates[56] - rev_rates[56]);
  //sp 11
  sp_rates[10] -= (fwd_rates[56] - rev_rates[56]);

  //rxn 57
  //sp 11
  sp_rates[10] -= (fwd_rates[57] - rev_rates[57]);
  //sp 10
  sp_rates[9] += (fwd_rates[57] - rev_rates[57]);

  //rxn 58
  //sp 17
  sp_rates[16] += (fwd_rates[58] - rev_rates[58]);
  //sp 18
  sp_rates[17] -= (fwd_rates[58] - rev_rates[58]);
  //sp 11
  sp_rates[10] -= (fwd_rates[58] - rev_rates[58]);
  //sp 15
  sp_rates[14] += (fwd_rates[58] - rev_rates[58]);

  //rxn 59
  //sp 2
  sp_rates[1] -= (fwd_rates[59] - rev_rates[59]) * pres_mod[11];
  //sp 15
  sp_rates[14] -= (fwd_rates[59] - rev_rates[59]) * pres_mod[11];
  //sp 16
  sp_rates[15] = (fwd_rates[59] - rev_rates[59]) * pres_mod[11];

  //rxn 60
  //sp 2
  sp_rates[1] -= (fwd_rates[60] - rev_rates[60]);
  //sp 6
  sp_rates[5] += (fwd_rates[60] - rev_rates[60]);
  //sp 15
  sp_rates[14] -= (fwd_rates[60] - rev_rates[60]);
  //sp 14
  sp_rates[13] += (fwd_rates[60] - rev_rates[60]);

  //rxn 61
  //sp 3
  sp_rates[2] -= (fwd_rates[61] - rev_rates[61]);
  //sp 4
  sp_rates[3] += (fwd_rates[61] - rev_rates[61]);
  //sp 14
  sp_rates[13] += (fwd_rates[61] - rev_rates[61]);
  //sp 15
  sp_rates[14] -= (fwd_rates[61] - rev_rates[61]);

  //rxn 62
  //sp 4
  sp_rates[3] -= (fwd_rates[62] - rev_rates[62]);
  //sp 7
  sp_rates[6] += (fwd_rates[62] - rev_rates[62]);
  //sp 15
  sp_rates[14] -= (fwd_rates[62] - rev_rates[62]);
  //sp 14
  sp_rates[13] += (fwd_rates[62] - rev_rates[62]);

  //rxn 63
  //sp 9
  sp_rates[8] -= (fwd_rates[63] - rev_rates[63]);
  //sp 5
  sp_rates[4] += (fwd_rates[63] - rev_rates[63]);
  //sp 14
  sp_rates[13] += (fwd_rates[63] - rev_rates[63]);
  //sp 15
  sp_rates[14] -= (fwd_rates[63] - rev_rates[63]);

  //rxn 64
  //sp 5
  sp_rates[4] -= (fwd_rates[64] - rev_rates[64]);
  //sp 14
  sp_rates[13] += (fwd_rates[64] - rev_rates[64]);
  //sp 15
  sp_rates[14] -= (fwd_rates[64] - rev_rates[64]);
  //sp 8
  sp_rates[7] += (fwd_rates[64] - rev_rates[64]);

  //rxn 65
  //sp 2
  sp_rates[1] -= (fwd_rates[65] - rev_rates[65]) * pres_mod[12];
  //sp 12
  sp_rates[11] -= (fwd_rates[65] - rev_rates[65]) * pres_mod[12];
  //sp 13
  sp_rates[12] = (fwd_rates[65] - rev_rates[65]) * pres_mod[12];

  //rxn 66
  //sp 2
  sp_rates[1] += (fwd_rates[66] - rev_rates[66]);
  //sp 3
  sp_rates[2] -= (fwd_rates[66] - rev_rates[66]);
  //sp 12
  sp_rates[11] -= (fwd_rates[66] - rev_rates[66]);
  //sp 15
  sp_rates[14] += (fwd_rates[66] - rev_rates[66]);

  //rxn 67
  //sp 4
  sp_rates[3] -= (fwd_rates[67] - rev_rates[67]);
  //sp 10
  sp_rates[9] += (fwd_rates[67] - rev_rates[67]);
  //sp 12
  sp_rates[11] -= (fwd_rates[67] - rev_rates[67]);
  //sp 7
  sp_rates[6] += (fwd_rates[67] - rev_rates[67]);

  //rxn 68
  //sp 4
  sp_rates[3] -= (fwd_rates[68] - rev_rates[68]);
  //sp 11
  sp_rates[10] += (fwd_rates[68] - rev_rates[68]);
  //sp 12
  sp_rates[11] -= (fwd_rates[68] - rev_rates[68]);
  //sp 7
  sp_rates[6] += (fwd_rates[68] - rev_rates[68]);

  //rxn 69
  //sp 9
  sp_rates[8] -= (fwd_rates[69] - rev_rates[69]);
  //sp 3
  sp_rates[2] += (fwd_rates[69] - rev_rates[69]);
  //sp 12
  sp_rates[11] -= (fwd_rates[69] - rev_rates[69]);
  //sp 16
  sp_rates[15] += (fwd_rates[69] - rev_rates[69]);

  //rxn 70
  //sp 9
  sp_rates[8] -= (fwd_rates[70] - rev_rates[70]);
  //sp 4
  sp_rates[3] += (fwd_rates[70] - rev_rates[70]);
  //sp 12
  sp_rates[11] -= (fwd_rates[70] - rev_rates[70]);
  //sp 15
  sp_rates[14] += (fwd_rates[70] - rev_rates[70]);

  //rxn 71
  //sp 9
  sp_rates[8] += (fwd_rates[71] - rev_rates[71]);
  //sp 13
  sp_rates[12] += (fwd_rates[71] - rev_rates[71]);
  //sp 12
  sp_rates[11] -= (fwd_rates[71] - rev_rates[71]);
  //sp 5
  sp_rates[4] -= (fwd_rates[71] - rev_rates[71]);

  //rxn 72
  //sp 4
  sp_rates[3] += (fwd_rates[72] - rev_rates[72]);
  //sp 12
  sp_rates[11] -= (fwd_rates[72] - rev_rates[72]);
  //sp 5
  sp_rates[4] -= (fwd_rates[72] - rev_rates[72]);
  //sp 16
  sp_rates[15] += (fwd_rates[72] - rev_rates[72]);

  //rxn 73
  //sp 5
  sp_rates[4] += (fwd_rates[73] - rev_rates[73]);
  //sp 12
  sp_rates[11] -= (fwd_rates[73] - rev_rates[73]);
  //sp 13
  sp_rates[12] += (fwd_rates[73] - rev_rates[73]);
  //sp 8
  sp_rates[7] -= (fwd_rates[73] - rev_rates[73]);

  //rxn 74
  //sp 17
  sp_rates[16] += (fwd_rates[74] - rev_rates[74]);
  //sp 12
  sp_rates[11] -= (fwd_rates[74] - rev_rates[74]);
  //sp 13
  sp_rates[12] += (fwd_rates[74] - rev_rates[74]);
  //sp 14
  sp_rates[13] -= (fwd_rates[74] - rev_rates[74]);

  //rxn 75
  //sp 12
  sp_rates[11] -= (fwd_rates[75] - rev_rates[75]);
  //sp 13
  sp_rates[12] += (fwd_rates[75] - rev_rates[75]);
  //sp 14
  sp_rates[13] += (fwd_rates[75] - rev_rates[75]);
  //sp 15
  sp_rates[14] -= (fwd_rates[75] - rev_rates[75]);

  //rxn 76
  //sp 10
  sp_rates[9] -= (fwd_rates[76] - rev_rates[76]);
  //sp 12
  sp_rates[11] -= (fwd_rates[76] - rev_rates[76]);
  //sp 21
  sp_rates[20] = (fwd_rates[76] - rev_rates[76]);
  //sp 2
  sp_rates[1] += (fwd_rates[76] - rev_rates[76]);

  //rxn 77
  //sp 2
  sp_rates[1] += (fwd_rates[77] - rev_rates[77]);
  //sp 11
  sp_rates[10] -= (fwd_rates[77] - rev_rates[77]);
  //sp 12
  sp_rates[11] -= (fwd_rates[77] - rev_rates[77]);
  //sp 21
  sp_rates[20] += (fwd_rates[77] - rev_rates[77]);

  //rxn 78
  //sp 12
  sp_rates[11] -= 2.0 * (fwd_rates[78] - rev_rates[78]) * pres_mod[13];
  //sp 23
  sp_rates[22] = (fwd_rates[78] - rev_rates[78]) * pres_mod[13];

  //rxn 79
  //sp 2
  sp_rates[1] += (fwd_rates[79] - rev_rates[79]);
  //sp 12
  sp_rates[11] -= 2.0 * (fwd_rates[79] - rev_rates[79]);
  //sp 22
  sp_rates[21] = (fwd_rates[79] - rev_rates[79]);

  //rxn 80
  //sp 2
  sp_rates[1] -= (fwd_rates[80] - rev_rates[80]);
  //sp 6
  sp_rates[5] += (fwd_rates[80] - rev_rates[80]);
  //sp 15
  sp_rates[14] += (fwd_rates[80] - rev_rates[80]);
  //sp 16
  sp_rates[15] -= (fwd_rates[80] - rev_rates[80]);

  //rxn 81
  //sp 4
  sp_rates[3] += (fwd_rates[81] - rev_rates[81]);
  //sp 2
  sp_rates[1] -= (fwd_rates[81] - rev_rates[81]);
  //sp 12
  sp_rates[11] += (fwd_rates[81] - rev_rates[81]);
  //sp 16
  sp_rates[15] -= (fwd_rates[81] - rev_rates[81]);

  //rxn 82
  //sp 2
  sp_rates[1] -= (fwd_rates[82] - rev_rates[82]);
  //sp 11
  sp_rates[10] += (fwd_rates[82] - rev_rates[82]);
  //sp 7
  sp_rates[6] += (fwd_rates[82] - rev_rates[82]);
  //sp 16
  sp_rates[15] -= (fwd_rates[82] - rev_rates[82]);

  //rxn 83
  //sp 3
  sp_rates[2] -= (fwd_rates[83] - rev_rates[83]);
  //sp 4
  sp_rates[3] += (fwd_rates[83] - rev_rates[83]);
  //sp 15
  sp_rates[14] += (fwd_rates[83] - rev_rates[83]);
  //sp 16
  sp_rates[15] -= (fwd_rates[83] - rev_rates[83]);

  //rxn 84
  //sp 4
  sp_rates[3] -= (fwd_rates[84] - rev_rates[84]);
  //sp 7
  sp_rates[6] += (fwd_rates[84] - rev_rates[84]);
  //sp 15
  sp_rates[14] += (fwd_rates[84] - rev_rates[84]);
  //sp 16
  sp_rates[15] -= (fwd_rates[84] - rev_rates[84]);

  //rxn 85
  //sp 9
  sp_rates[8] -= (fwd_rates[85] - rev_rates[85]);
  //sp 5
  sp_rates[4] += (fwd_rates[85] - rev_rates[85]);
  //sp 15
  sp_rates[14] += (fwd_rates[85] - rev_rates[85]);
  //sp 16
  sp_rates[15] -= (fwd_rates[85] - rev_rates[85]);

  //rxn 86
  //sp 2
  sp_rates[1] -= (fwd_rates[86] - rev_rates[86]);
  //sp 12
  sp_rates[11] += (fwd_rates[86] - rev_rates[86]);
  //sp 13
  sp_rates[12] -= (fwd_rates[86] - rev_rates[86]);
  //sp 6
  sp_rates[5] += (fwd_rates[86] - rev_rates[86]);

  //rxn 87
  //sp 4
  sp_rates[3] += (fwd_rates[87] - rev_rates[87]);
  //sp 3
  sp_rates[2] -= (fwd_rates[87] - rev_rates[87]);
  //sp 12
  sp_rates[11] += (fwd_rates[87] - rev_rates[87]);
  //sp 13
  sp_rates[12] -= (fwd_rates[87] - rev_rates[87]);

  //rxn 88
  //sp 12
  sp_rates[11] += (fwd_rates[88] - rev_rates[88]);
  //sp 4
  sp_rates[3] -= (fwd_rates[88] - rev_rates[88]);
  //sp 13
  sp_rates[12] -= (fwd_rates[88] - rev_rates[88]);
  //sp 7
  sp_rates[6] += (fwd_rates[88] - rev_rates[88]);

  //rxn 89
  //sp 10
  sp_rates[9] -= (fwd_rates[89] - rev_rates[89]);
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[89] - rev_rates[89]);
  //sp 13
  sp_rates[12] -= (fwd_rates[89] - rev_rates[89]);

  //rxn 90
  //sp 11
  sp_rates[10] -= (fwd_rates[90] - rev_rates[90]);
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[90] - rev_rates[90]);
  //sp 13
  sp_rates[12] -= (fwd_rates[90] - rev_rates[90]);

  //rxn 91
  //sp 2
  sp_rates[1] += (fwd_rates[91] - rev_rates[91]) * pres_mod[14];
  //sp 19
  sp_rates[18] += (fwd_rates[91] - rev_rates[91]) * pres_mod[14];
  //sp 20
  sp_rates[19] = -(fwd_rates[91] - rev_rates[91]) * pres_mod[14];

  //rxn 92
  //sp 17
  sp_rates[16] += (fwd_rates[92] - rev_rates[92]);
  //sp 10
  sp_rates[9] += (fwd_rates[92] - rev_rates[92]);
  //sp 19
  sp_rates[18] -= (fwd_rates[92] - rev_rates[92]);
  //sp 3
  sp_rates[2] -= (fwd_rates[92] - rev_rates[92]);

  //rxn 93
  //sp 12
  sp_rates[11] += (fwd_rates[93] - rev_rates[93]);
  //sp 17
  sp_rates[16] += (fwd_rates[93] - rev_rates[93]);
  //sp 19
  sp_rates[18] -= (fwd_rates[93] - rev_rates[93]);
  //sp 4
  sp_rates[3] -= (fwd_rates[93] - rev_rates[93]);

  //rxn 94
  //sp 17
  sp_rates[16] += (fwd_rates[94] - rev_rates[94]);
  //sp 19
  sp_rates[18] -= (fwd_rates[94] - rev_rates[94]);
  //sp 20
  sp_rates[19] += (fwd_rates[94] - rev_rates[94]);
  //sp 14
  sp_rates[13] -= (fwd_rates[94] - rev_rates[94]);

  //rxn 95
  //sp 25
  sp_rates[24] = (fwd_rates[95] - rev_rates[95]);
  //sp 19
  sp_rates[18] -= (fwd_rates[95] - rev_rates[95]);
  //sp 12
  sp_rates[11] -= (fwd_rates[95] - rev_rates[95]);

  //rxn 96
  //sp 2
  sp_rates[1] -= (fwd_rates[96] - rev_rates[96]) * pres_mod[15];
  //sp 20
  sp_rates[19] -= (fwd_rates[96] - rev_rates[96]) * pres_mod[15];
  //sp 21
  sp_rates[20] += (fwd_rates[96] - rev_rates[96]) * pres_mod[15];

  //rxn 97
  //sp 2
  sp_rates[1] -= (fwd_rates[97] - rev_rates[97]);
  //sp 19
  sp_rates[18] += (fwd_rates[97] - rev_rates[97]);
  //sp 20
  sp_rates[19] -= (fwd_rates[97] - rev_rates[97]);
  //sp 6
  sp_rates[5] += (fwd_rates[97] - rev_rates[97]);

  //rxn 98
  //sp 12
  sp_rates[11] += (fwd_rates[98] - rev_rates[98]);
  //sp 17
  sp_rates[16] += (fwd_rates[98] - rev_rates[98]);
  //sp 3
  sp_rates[2] -= (fwd_rates[98] - rev_rates[98]);
  //sp 20
  sp_rates[19] -= (fwd_rates[98] - rev_rates[98]);

  //rxn 99
  //sp 4
  sp_rates[3] -= (fwd_rates[99] - rev_rates[99]);
  //sp 19
  sp_rates[18] += (fwd_rates[99] - rev_rates[99]);
  //sp 20
  sp_rates[19] -= (fwd_rates[99] - rev_rates[99]);
  //sp 7
  sp_rates[6] += (fwd_rates[99] - rev_rates[99]);

  //rxn 100
  //sp 9
  sp_rates[8] -= (fwd_rates[100] - rev_rates[100]);
  //sp 19
  sp_rates[18] += (fwd_rates[100] - rev_rates[100]);
  //sp 20
  sp_rates[19] -= (fwd_rates[100] - rev_rates[100]);
  //sp 5
  sp_rates[4] += (fwd_rates[100] - rev_rates[100]);

  //rxn 101
  //sp 9
  sp_rates[8] -= (fwd_rates[101] - rev_rates[101]);
  //sp 3
  sp_rates[2] += (fwd_rates[101] - rev_rates[101]);
  //sp 20
  sp_rates[19] -= (fwd_rates[101] - rev_rates[101]);
  //sp 24
  sp_rates[23] = (fwd_rates[101] - rev_rates[101]);

  //rxn 102
  //sp 9
  sp_rates[8] -= (fwd_rates[102] - rev_rates[102]);
  //sp 20
  sp_rates[19] -= (fwd_rates[102] - rev_rates[102]);
  //sp 14
  sp_rates[13] += (fwd_rates[102] - rev_rates[102]);
  //sp 15
  sp_rates[14] += (fwd_rates[102] - rev_rates[102]);

  //rxn 103
  //sp 4
  sp_rates[3] += (fwd_rates[103] - rev_rates[103]);
  //sp 20
  sp_rates[19] -= (fwd_rates[103] - rev_rates[103]);
  //sp 5
  sp_rates[4] -= (fwd_rates[103] - rev_rates[103]);
  //sp 24
  sp_rates[23] += (fwd_rates[103] - rev_rates[103]);

  //rxn 104
  //sp 5
  sp_rates[4] += (fwd_rates[104] - rev_rates[104]);
  //sp 20
  sp_rates[19] -= (fwd_rates[104] - rev_rates[104]);
  //sp 21
  sp_rates[20] += (fwd_rates[104] - rev_rates[104]);
  //sp 8
  sp_rates[7] -= (fwd_rates[104] - rev_rates[104]);

  //rxn 105
  //sp 17
  sp_rates[16] += (fwd_rates[105] - rev_rates[105]);
  //sp 20
  sp_rates[19] -= (fwd_rates[105] - rev_rates[105]);
  //sp 21
  sp_rates[20] += (fwd_rates[105] - rev_rates[105]);
  //sp 14
  sp_rates[13] -= (fwd_rates[105] - rev_rates[105]);

  //rxn 106
  //sp 28
  sp_rates[27] = (fwd_rates[106] - rev_rates[106]);
  //sp 20
  sp_rates[19] -= (fwd_rates[106] - rev_rates[106]);
  //sp 14
  sp_rates[13] -= (fwd_rates[106] - rev_rates[106]);

  //rxn 107
  //sp 20
  sp_rates[19] -= (fwd_rates[107] - rev_rates[107]);
  //sp 19
  sp_rates[18] += (fwd_rates[107] - rev_rates[107]);
  //sp 12
  sp_rates[11] -= (fwd_rates[107] - rev_rates[107]);
  //sp 13
  sp_rates[12] += (fwd_rates[107] - rev_rates[107]);

  //rxn 108
  //sp 20
  sp_rates[19] -= (fwd_rates[108] - rev_rates[108]) * pres_mod[16];
  //sp 26
  sp_rates[25] = (fwd_rates[108] - rev_rates[108]) * pres_mod[16];
  //sp 12
  sp_rates[11] -= (fwd_rates[108] - rev_rates[108]) * pres_mod[16];

  //rxn 109
  //sp 20
  sp_rates[19] -= (fwd_rates[109] - rev_rates[109]);
  //sp 2
  sp_rates[1] += (fwd_rates[109] - rev_rates[109]);
  //sp 12
  sp_rates[11] -= (fwd_rates[109] - rev_rates[109]);
  //sp 25
  sp_rates[24] += (fwd_rates[109] - rev_rates[109]);

  //rxn 110
  //sp 19
  sp_rates[18] += (fwd_rates[110] - rev_rates[110]);
  //sp 20
  sp_rates[19] -= 2.0 * (fwd_rates[110] - rev_rates[110]);
  //sp 21
  sp_rates[20] += (fwd_rates[110] - rev_rates[110]);

  //rxn 111
  //sp 17
  sp_rates[16] += (fwd_rates[111] - rev_rates[111]);
  //sp 12
  sp_rates[11] += (fwd_rates[111] - rev_rates[111]);
  //sp 24
  sp_rates[23] -= (fwd_rates[111] - rev_rates[111]);

  //rxn 112
  //sp 2
  sp_rates[1] -= (fwd_rates[112] - rev_rates[112]);
  //sp 12
  sp_rates[11] += (fwd_rates[112] - rev_rates[112]);
  //sp 14
  sp_rates[13] += (fwd_rates[112] - rev_rates[112]);
  //sp 24
  sp_rates[23] -= (fwd_rates[112] - rev_rates[112]);

  //rxn 113
  //sp 9
  sp_rates[8] -= (fwd_rates[113] - rev_rates[113]);
  //sp 17
  sp_rates[16] += (fwd_rates[113] - rev_rates[113]);
  //sp 4
  sp_rates[3] += (fwd_rates[113] - rev_rates[113]);
  //sp 15
  sp_rates[14] += (fwd_rates[113] - rev_rates[113]);
  //sp 24
  sp_rates[23] -= (fwd_rates[113] - rev_rates[113]);

  //rxn 114
  //sp 2
  sp_rates[1] -= (fwd_rates[114] - rev_rates[114]) * pres_mod[17];
  //sp 21
  sp_rates[20] -= (fwd_rates[114] - rev_rates[114]) * pres_mod[17];
  //sp 22
  sp_rates[21] += (fwd_rates[114] - rev_rates[114]) * pres_mod[17];

  //rxn 115
  //sp 2
  sp_rates[1] -= (fwd_rates[115] - rev_rates[115]);
  //sp 20
  sp_rates[19] += (fwd_rates[115] - rev_rates[115]);
  //sp 21
  sp_rates[20] -= (fwd_rates[115] - rev_rates[115]);
  //sp 6
  sp_rates[5] += (fwd_rates[115] - rev_rates[115]);

  //rxn 116
  //sp 4
  sp_rates[3] += (fwd_rates[116] - rev_rates[116]);
  //sp 3
  sp_rates[2] -= (fwd_rates[116] - rev_rates[116]);
  //sp 20
  sp_rates[19] += (fwd_rates[116] - rev_rates[116]);
  //sp 21
  sp_rates[20] -= (fwd_rates[116] - rev_rates[116]);

  //rxn 117
  //sp 3
  sp_rates[2] -= (fwd_rates[117] - rev_rates[117]);
  //sp 12
  sp_rates[11] += (fwd_rates[117] - rev_rates[117]);
  //sp 21
  sp_rates[20] -= (fwd_rates[117] - rev_rates[117]);
  //sp 14
  sp_rates[13] += (fwd_rates[117] - rev_rates[117]);

  //rxn 118
  //sp 10
  sp_rates[9] += (fwd_rates[118] - rev_rates[118]);
  //sp 3
  sp_rates[2] -= (fwd_rates[118] - rev_rates[118]);
  //sp 21
  sp_rates[20] -= (fwd_rates[118] - rev_rates[118]);
  //sp 15
  sp_rates[14] += (fwd_rates[118] - rev_rates[118]);

  //rxn 119
  //sp 20
  sp_rates[19] += (fwd_rates[119] - rev_rates[119]);
  //sp 4
  sp_rates[3] -= (fwd_rates[119] - rev_rates[119]);
  //sp 21
  sp_rates[20] -= (fwd_rates[119] - rev_rates[119]);
  //sp 7
  sp_rates[6] += (fwd_rates[119] - rev_rates[119]);

  //rxn 120
  //sp 17
  sp_rates[16] += (fwd_rates[120] - rev_rates[120]);
  //sp 21
  sp_rates[20] -= (fwd_rates[120] - rev_rates[120]);
  //sp 14
  sp_rates[13] -= (fwd_rates[120] - rev_rates[120]);
  //sp 22
  sp_rates[21] += (fwd_rates[120] - rev_rates[120]);

  //rxn 121
  //sp 25
  sp_rates[24] += (fwd_rates[121] - rev_rates[121]);
  //sp 10
  sp_rates[9] -= (fwd_rates[121] - rev_rates[121]);
  //sp 21
  sp_rates[20] -= (fwd_rates[121] - rev_rates[121]);
  //sp 2
  sp_rates[1] += (fwd_rates[121] - rev_rates[121]);

  //rxn 122
  //sp 25
  sp_rates[24] += (fwd_rates[122] - rev_rates[122]);
  //sp 2
  sp_rates[1] += (fwd_rates[122] - rev_rates[122]);
  //sp 11
  sp_rates[10] -= (fwd_rates[122] - rev_rates[122]);
  //sp 21
  sp_rates[20] -= (fwd_rates[122] - rev_rates[122]);

  //rxn 123
  //sp 20
  sp_rates[19] += (fwd_rates[123] - rev_rates[123]);
  //sp 13
  sp_rates[12] += (fwd_rates[123] - rev_rates[123]);
  //sp 12
  sp_rates[11] -= (fwd_rates[123] - rev_rates[123]);
  //sp 21
  sp_rates[20] -= (fwd_rates[123] - rev_rates[123]);

  //rxn 124
  //sp 27
  sp_rates[26] = -(fwd_rates[124] - rev_rates[124]);
  //sp 12
  sp_rates[11] += (fwd_rates[124] - rev_rates[124]);
  //sp 21
  sp_rates[20] += (fwd_rates[124] - rev_rates[124]);

  //rxn 125
  //sp 9
  sp_rates[8] -= (fwd_rates[125] - rev_rates[125]);
  //sp 5
  sp_rates[4] += (fwd_rates[125] - rev_rates[125]);
  //sp 20
  sp_rates[19] += (fwd_rates[125] - rev_rates[125]);
  //sp 21
  sp_rates[20] -= (fwd_rates[125] - rev_rates[125]);

  //rxn 126
  //sp 29
  sp_rates[28] = (fwd_rates[126] - rev_rates[126]);
  //sp 20
  sp_rates[19] -= (fwd_rates[126] - rev_rates[126]);
  //sp 21
  sp_rates[20] -= (fwd_rates[126] - rev_rates[126]);

  //rxn 127
  //sp 2
  sp_rates[1] -= (fwd_rates[127] - rev_rates[127]) * pres_mod[18];
  //sp 22
  sp_rates[21] -= (fwd_rates[127] - rev_rates[127]) * pres_mod[18];
  //sp 23
  sp_rates[22] += (fwd_rates[127] - rev_rates[127]) * pres_mod[18];

  //rxn 128
  //sp 2
  sp_rates[1] -= (fwd_rates[128] - rev_rates[128]);
  //sp 21
  sp_rates[20] += (fwd_rates[128] - rev_rates[128]);
  //sp 22
  sp_rates[21] -= (fwd_rates[128] - rev_rates[128]);
  //sp 6
  sp_rates[5] += (fwd_rates[128] - rev_rates[128]);

  //rxn 129
  //sp 3
  sp_rates[2] -= (fwd_rates[129] - rev_rates[129]);
  //sp 12
  sp_rates[11] += (fwd_rates[129] - rev_rates[129]);
  //sp 22
  sp_rates[21] -= (fwd_rates[129] - rev_rates[129]);
  //sp 15
  sp_rates[14] += (fwd_rates[129] - rev_rates[129]);

  //rxn 130
  //sp 9
  sp_rates[8] -= (fwd_rates[130] - rev_rates[130]);
  //sp 5
  sp_rates[4] += (fwd_rates[130] - rev_rates[130]);
  //sp 21
  sp_rates[20] += (fwd_rates[130] - rev_rates[130]);
  //sp 22
  sp_rates[21] -= (fwd_rates[130] - rev_rates[130]);

  //rxn 131
  //sp 9
  sp_rates[8] += (fwd_rates[131] - rev_rates[131]);
  //sp 5
  sp_rates[4] -= (fwd_rates[131] - rev_rates[131]);
  //sp 22
  sp_rates[21] -= (fwd_rates[131] - rev_rates[131]);
  //sp 23
  sp_rates[22] += (fwd_rates[131] - rev_rates[131]);

  //rxn 132
  //sp 21
  sp_rates[20] += (fwd_rates[132] - rev_rates[132]);
  //sp 5
  sp_rates[4] -= (fwd_rates[132] - rev_rates[132]);
  //sp 22
  sp_rates[21] -= (fwd_rates[132] - rev_rates[132]);
  //sp 8
  sp_rates[7] += (fwd_rates[132] - rev_rates[132]);

  //rxn 133
  //sp 4
  sp_rates[3] += (fwd_rates[133] - rev_rates[133]);
  //sp 12
  sp_rates[11] += (fwd_rates[133] - rev_rates[133]);
  //sp 5
  sp_rates[4] -= (fwd_rates[133] - rev_rates[133]);
  //sp 22
  sp_rates[21] -= (fwd_rates[133] - rev_rates[133]);
  //sp 15
  sp_rates[14] += (fwd_rates[133] - rev_rates[133]);

  //rxn 134
  //sp 5
  sp_rates[4] += (fwd_rates[134] - rev_rates[134]);
  //sp 22
  sp_rates[21] -= (fwd_rates[134] - rev_rates[134]);
  //sp 23
  sp_rates[22] += (fwd_rates[134] - rev_rates[134]);
  //sp 8
  sp_rates[7] -= (fwd_rates[134] - rev_rates[134]);

  //rxn 135
  //sp 20
  sp_rates[19] -= (fwd_rates[135] - rev_rates[135]) * pres_mod[19];
  //sp 22
  sp_rates[21] -= (fwd_rates[135] - rev_rates[135]) * pres_mod[19];
  //sp 30
  sp_rates[29] = (fwd_rates[135] - rev_rates[135]) * pres_mod[19];

  //rxn 136
  //sp 12
  sp_rates[11] += (fwd_rates[136] - rev_rates[136]);
  //sp 25
  sp_rates[24] += (fwd_rates[136] - rev_rates[136]);
  //sp 20
  sp_rates[19] -= (fwd_rates[136] - rev_rates[136]);
  //sp 22
  sp_rates[21] -= (fwd_rates[136] - rev_rates[136]);

  //rxn 137
  //sp 2
  sp_rates[1] -= (fwd_rates[137] - rev_rates[137]);
  //sp 6
  sp_rates[5] += (fwd_rates[137] - rev_rates[137]);
  //sp 23
  sp_rates[22] -= (fwd_rates[137] - rev_rates[137]);
  //sp 22
  sp_rates[21] += (fwd_rates[137] - rev_rates[137]);

  //rxn 138
  //sp 3
  sp_rates[2] -= (fwd_rates[138] - rev_rates[138]);
  //sp 4
  sp_rates[3] += (fwd_rates[138] - rev_rates[138]);
  //sp 22
  sp_rates[21] += (fwd_rates[138] - rev_rates[138]);
  //sp 23
  sp_rates[22] -= (fwd_rates[138] - rev_rates[138]);

  //rxn 139
  //sp 7
  sp_rates[6] += (fwd_rates[139] - rev_rates[139]);
  //sp 4
  sp_rates[3] -= (fwd_rates[139] - rev_rates[139]);
  //sp 22
  sp_rates[21] += (fwd_rates[139] - rev_rates[139]);
  //sp 23
  sp_rates[22] -= (fwd_rates[139] - rev_rates[139]);

  //rxn 140
  //sp 11
  sp_rates[10] -= (fwd_rates[140] - rev_rates[140]);
  //sp 12
  sp_rates[11] += (fwd_rates[140] - rev_rates[140]);
  //sp 22
  sp_rates[21] += (fwd_rates[140] - rev_rates[140]);
  //sp 23
  sp_rates[22] -= (fwd_rates[140] - rev_rates[140]);

  //rxn 141
  //sp 12
  sp_rates[11] -= (fwd_rates[141] - rev_rates[141]);
  //sp 13
  sp_rates[12] += (fwd_rates[141] - rev_rates[141]);
  //sp 22
  sp_rates[21] += (fwd_rates[141] - rev_rates[141]);
  //sp 23
  sp_rates[22] -= (fwd_rates[141] - rev_rates[141]);

  //rxn 142
  //sp 25
  sp_rates[24] -= (fwd_rates[142] - rev_rates[142]) * pres_mod[20];
  //sp 2
  sp_rates[1] -= (fwd_rates[142] - rev_rates[142]) * pres_mod[20];
  //sp 26
  sp_rates[25] += (fwd_rates[142] - rev_rates[142]) * pres_mod[20];

  //rxn 143
  //sp 25
  sp_rates[24] -= (fwd_rates[143] - rev_rates[143]);
  //sp 2
  sp_rates[1] += (fwd_rates[143] - rev_rates[143]);
  //sp 3
  sp_rates[2] -= (fwd_rates[143] - rev_rates[143]);
  //sp 28
  sp_rates[27] += (fwd_rates[143] - rev_rates[143]);

  //rxn 144
  //sp 25
  sp_rates[24] -= (fwd_rates[144] - rev_rates[144]);
  //sp 2
  sp_rates[1] += 2.0 * (fwd_rates[144] - rev_rates[144]);
  //sp 4
  sp_rates[3] -= (fwd_rates[144] - rev_rates[144]);
  //sp 28
  sp_rates[27] += (fwd_rates[144] - rev_rates[144]);

  //rxn 145
  //sp 25
  sp_rates[24] -= (fwd_rates[145] - rev_rates[145]);
  //sp 9
  sp_rates[8] -= (fwd_rates[145] - rev_rates[145]);
  //sp 28
  sp_rates[27] += (fwd_rates[145] - rev_rates[145]);
  //sp 4
  sp_rates[3] += (fwd_rates[145] - rev_rates[145]);

  //rxn 146
  //sp 25
  sp_rates[24] -= (fwd_rates[146] - rev_rates[146]);
  //sp 26
  sp_rates[25] += (fwd_rates[146] - rev_rates[146]);
  //sp 5
  sp_rates[4] -= (fwd_rates[146] - rev_rates[146]);
  //sp 9
  sp_rates[8] += (fwd_rates[146] - rev_rates[146]);

  //rxn 147
  //sp 25
  sp_rates[24] -= (fwd_rates[147] - rev_rates[147]);
  //sp 4
  sp_rates[3] += (fwd_rates[147] - rev_rates[147]);
  //sp 20
  sp_rates[19] += (fwd_rates[147] - rev_rates[147]);
  //sp 5
  sp_rates[4] -= (fwd_rates[147] - rev_rates[147]);
  //sp 15
  sp_rates[14] += (fwd_rates[147] - rev_rates[147]);

  //rxn 148
  //sp 25
  sp_rates[24] -= (fwd_rates[148] - rev_rates[148]);
  //sp 26
  sp_rates[25] += (fwd_rates[148] - rev_rates[148]);
  //sp 14
  sp_rates[13] -= (fwd_rates[148] - rev_rates[148]);
  //sp 17
  sp_rates[16] += (fwd_rates[148] - rev_rates[148]);

  //rxn 149
  //sp 25
  sp_rates[24] -= (fwd_rates[149] - rev_rates[149]) * pres_mod[21];
  //sp 12
  sp_rates[11] -= (fwd_rates[149] - rev_rates[149]) * pres_mod[21];
  //sp 30
  sp_rates[29] += (fwd_rates[149] - rev_rates[149]) * pres_mod[21];

  //rxn 150
  //sp 26
  sp_rates[25] -= (fwd_rates[150] - rev_rates[150]) * pres_mod[22];
  //sp 27
  sp_rates[26] += (fwd_rates[150] - rev_rates[150]) * pres_mod[22];
  //sp 2
  sp_rates[1] -= (fwd_rates[150] - rev_rates[150]) * pres_mod[22];

  //rxn 151
  //sp 26
  sp_rates[25] -= (fwd_rates[151] - rev_rates[151]);
  //sp 12
  sp_rates[11] += (fwd_rates[151] - rev_rates[151]);
  //sp 21
  sp_rates[20] += (fwd_rates[151] - rev_rates[151]);
  //sp 2
  sp_rates[1] -= (fwd_rates[151] - rev_rates[151]);

  //rxn 152
  //sp 25
  sp_rates[24] += (fwd_rates[152] - rev_rates[152]);
  //sp 26
  sp_rates[25] -= (fwd_rates[152] - rev_rates[152]);
  //sp 6
  sp_rates[5] += (fwd_rates[152] - rev_rates[152]);
  //sp 2
  sp_rates[1] -= (fwd_rates[152] - rev_rates[152]);

  //rxn 153
  //sp 26
  sp_rates[25] -= (fwd_rates[153] - rev_rates[153]);
  //sp 3
  sp_rates[2] -= (fwd_rates[153] - rev_rates[153]);
  //sp 28
  sp_rates[27] += (fwd_rates[153] - rev_rates[153]);
  //sp 2
  sp_rates[1] += 2.0 * (fwd_rates[153] - rev_rates[153]);

  //rxn 154
  //sp 26
  sp_rates[25] -= (fwd_rates[154] - rev_rates[154]);
  //sp 3
  sp_rates[2] -= (fwd_rates[154] - rev_rates[154]);
  //sp 22
  sp_rates[21] += (fwd_rates[154] - rev_rates[154]);
  //sp 14
  sp_rates[13] += (fwd_rates[154] - rev_rates[154]);

  //rxn 155
  //sp 25
  sp_rates[24] += (fwd_rates[155] - rev_rates[155]);
  //sp 26
  sp_rates[25] -= (fwd_rates[155] - rev_rates[155]);
  //sp 3
  sp_rates[2] -= (fwd_rates[155] - rev_rates[155]);
  //sp 4
  sp_rates[3] += (fwd_rates[155] - rev_rates[155]);

  //rxn 156
  //sp 25
  sp_rates[24] += (fwd_rates[156] - rev_rates[156]);
  //sp 26
  sp_rates[25] -= (fwd_rates[156] - rev_rates[156]);
  //sp 4
  sp_rates[3] -= (fwd_rates[156] - rev_rates[156]);
  //sp 7
  sp_rates[6] += (fwd_rates[156] - rev_rates[156]);

  //rxn 157
  //sp 25
  sp_rates[24] += (fwd_rates[157] - rev_rates[157]);
  //sp 26
  sp_rates[25] -= (fwd_rates[157] - rev_rates[157]);
  //sp 5
  sp_rates[4] -= (fwd_rates[157] - rev_rates[157]);
  //sp 8
  sp_rates[7] += (fwd_rates[157] - rev_rates[157]);

  //rxn 158
  //sp 25
  sp_rates[24] += (fwd_rates[158] - rev_rates[158]);
  //sp 26
  sp_rates[25] -= (fwd_rates[158] - rev_rates[158]);
  //sp 12
  sp_rates[11] -= (fwd_rates[158] - rev_rates[158]);
  //sp 13
  sp_rates[12] += (fwd_rates[158] - rev_rates[158]);

  //rxn 159
  //sp 2
  sp_rates[1] -= (fwd_rates[159] - rev_rates[159]);
  //sp 28
  sp_rates[27] -= (fwd_rates[159] - rev_rates[159]);
  //sp 21
  sp_rates[20] += (fwd_rates[159] - rev_rates[159]);
  //sp 14
  sp_rates[13] += (fwd_rates[159] - rev_rates[159]);

  //rxn 160
  //sp 20
  sp_rates[19] += (fwd_rates[160] - rev_rates[160]);
  //sp 17
  sp_rates[16] += (fwd_rates[160] - rev_rates[160]);
  //sp 3
  sp_rates[2] -= (fwd_rates[160] - rev_rates[160]);
  //sp 28
  sp_rates[27] -= (fwd_rates[160] - rev_rates[160]);
  //sp 4
  sp_rates[3] += (fwd_rates[160] - rev_rates[160]);

  //rxn 161
  //sp 4
  sp_rates[3] -= (fwd_rates[161] - rev_rates[161]);
  //sp 20
  sp_rates[19] += (fwd_rates[161] - rev_rates[161]);
  //sp 28
  sp_rates[27] -= (fwd_rates[161] - rev_rates[161]);
  //sp 7
  sp_rates[6] += (fwd_rates[161] - rev_rates[161]);
  //sp 17
  sp_rates[16] += (fwd_rates[161] - rev_rates[161]);

  //rxn 162
  //sp 2
  sp_rates[1] -= (fwd_rates[162] - rev_rates[162]);
  //sp 27
  sp_rates[26] -= (fwd_rates[162] - rev_rates[162]);
  //sp 12
  sp_rates[11] += (fwd_rates[162] - rev_rates[162]);
  //sp 22
  sp_rates[21] += (fwd_rates[162] - rev_rates[162]);

  //rxn 163
  //sp 2
  sp_rates[1] -= (fwd_rates[163] - rev_rates[163]);
  //sp 27
  sp_rates[26] -= (fwd_rates[163] - rev_rates[163]);
  //sp 6
  sp_rates[5] += (fwd_rates[163] - rev_rates[163]);
  //sp 26
  sp_rates[25] += (fwd_rates[163] - rev_rates[163]);

  //rxn 164
  //sp 27
  sp_rates[26] -= (fwd_rates[164] - rev_rates[164]);
  //sp 3
  sp_rates[2] -= (fwd_rates[164] - rev_rates[164]);
  //sp 15
  sp_rates[14] += (fwd_rates[164] - rev_rates[164]);
  //sp 22
  sp_rates[21] += (fwd_rates[164] - rev_rates[164]);

  //rxn 165
  //sp 26
  sp_rates[25] += (fwd_rates[165] - rev_rates[165]);
  //sp 27
  sp_rates[26] -= (fwd_rates[165] - rev_rates[165]);
  //sp 4
  sp_rates[3] -= (fwd_rates[165] - rev_rates[165]);
  //sp 7
  sp_rates[6] += (fwd_rates[165] - rev_rates[165]);

  //rxn 166
  //sp 9
  sp_rates[8] -= (fwd_rates[166] - rev_rates[166]);
  //sp 26
  sp_rates[25] += (fwd_rates[166] - rev_rates[166]);
  //sp 27
  sp_rates[26] -= (fwd_rates[166] - rev_rates[166]);
  //sp 5
  sp_rates[4] += (fwd_rates[166] - rev_rates[166]);

  //rxn 167
  //sp 27
  sp_rates[26] -= (fwd_rates[167] - rev_rates[167]);
  //sp 4
  sp_rates[3] += (fwd_rates[167] - rev_rates[167]);
  //sp 5
  sp_rates[4] -= (fwd_rates[167] - rev_rates[167]);
  //sp 22
  sp_rates[21] += (fwd_rates[167] - rev_rates[167]);
  //sp 15
  sp_rates[14] += (fwd_rates[167] - rev_rates[167]);

  //rxn 168
  //sp 26
  sp_rates[25] += (fwd_rates[168] - rev_rates[168]);
  //sp 27
  sp_rates[26] -= (fwd_rates[168] - rev_rates[168]);
  //sp 12
  sp_rates[11] -= (fwd_rates[168] - rev_rates[168]);
  //sp 13
  sp_rates[12] += (fwd_rates[168] - rev_rates[168]);

  //rxn 169
  //sp 2
  sp_rates[1] -= (fwd_rates[169] - rev_rates[169]) * pres_mod[23];
  //sp 29
  sp_rates[28] -= (fwd_rates[169] - rev_rates[169]) * pres_mod[23];
  //sp 30
  sp_rates[29] += (fwd_rates[169] - rev_rates[169]) * pres_mod[23];

  //rxn 170
  //sp 25
  sp_rates[24] += (fwd_rates[170] - rev_rates[170]);
  //sp 2
  sp_rates[1] -= (fwd_rates[170] - rev_rates[170]);
  //sp 12
  sp_rates[11] += (fwd_rates[170] - rev_rates[170]);
  //sp 29
  sp_rates[28] -= (fwd_rates[170] - rev_rates[170]);

  //rxn 171
  //sp 25
  sp_rates[24] += (fwd_rates[171] - rev_rates[171]);
  //sp 29
  sp_rates[28] -= (fwd_rates[171] - rev_rates[171]);
  //sp 4
  sp_rates[3] += (fwd_rates[171] - rev_rates[171]);
  //sp 5
  sp_rates[4] -= (fwd_rates[171] - rev_rates[171]);
  //sp 15
  sp_rates[14] += (fwd_rates[171] - rev_rates[171]);

  //rxn 172
  //sp 17
  sp_rates[16] += (fwd_rates[172] - rev_rates[172]);
  //sp 29
  sp_rates[28] -= (fwd_rates[172] - rev_rates[172]);
  //sp 14
  sp_rates[13] -= (fwd_rates[172] - rev_rates[172]);
  //sp 30
  sp_rates[29] += (fwd_rates[172] - rev_rates[172]);

  //rxn 173
  //sp 2
  sp_rates[1] -= (fwd_rates[173] - rev_rates[173]) * pres_mod[24];
  //sp 30
  sp_rates[29] -= (fwd_rates[173] - rev_rates[173]) * pres_mod[24];
  //sp 31
  sp_rates[30] = (fwd_rates[173] - rev_rates[173]) * pres_mod[24];

  //rxn 174
  //sp 2
  sp_rates[1] -= (fwd_rates[174] - rev_rates[174]);
  //sp 21
  sp_rates[20] += (fwd_rates[174] - rev_rates[174]);
  //sp 30
  sp_rates[29] -= (fwd_rates[174] - rev_rates[174]);
  //sp 22
  sp_rates[21] += (fwd_rates[174] - rev_rates[174]);

  //rxn 175
  //sp 2
  sp_rates[1] -= (fwd_rates[175] - rev_rates[175]);
  //sp 12
  sp_rates[11] += (fwd_rates[175] - rev_rates[175]);
  //sp 30
  sp_rates[29] -= (fwd_rates[175] - rev_rates[175]);
  //sp 26
  sp_rates[25] += (fwd_rates[175] - rev_rates[175]);

  //rxn 176
  //sp 2
  sp_rates[1] -= (fwd_rates[176] - rev_rates[176]);
  //sp 29
  sp_rates[28] += (fwd_rates[176] - rev_rates[176]);
  //sp 30
  sp_rates[29] -= (fwd_rates[176] - rev_rates[176]);
  //sp 6
  sp_rates[5] += (fwd_rates[176] - rev_rates[176]);

  //rxn 177
  //sp 3
  sp_rates[2] -= (fwd_rates[177] - rev_rates[177]);
  //sp 27
  sp_rates[26] += (fwd_rates[177] - rev_rates[177]);
  //sp 30
  sp_rates[29] -= (fwd_rates[177] - rev_rates[177]);
  //sp 14
  sp_rates[13] += (fwd_rates[177] - rev_rates[177]);

  //rxn 178
  //sp 3
  sp_rates[2] -= (fwd_rates[178] - rev_rates[178]);
  //sp 4
  sp_rates[3] += (fwd_rates[178] - rev_rates[178]);
  //sp 29
  sp_rates[28] += (fwd_rates[178] - rev_rates[178]);
  //sp 30
  sp_rates[29] -= (fwd_rates[178] - rev_rates[178]);

  //rxn 179
  //sp 3
  sp_rates[2] -= (fwd_rates[179] - rev_rates[179]);
  //sp 4
  sp_rates[3] += (fwd_rates[179] - rev_rates[179]);
  //sp 29
  sp_rates[28] += (fwd_rates[179] - rev_rates[179]);
  //sp 30
  sp_rates[29] -= (fwd_rates[179] - rev_rates[179]);

  //rxn 180
  //sp 4
  sp_rates[3] -= (fwd_rates[180] - rev_rates[180]);
  //sp 29
  sp_rates[28] += (fwd_rates[180] - rev_rates[180]);
  //sp 30
  sp_rates[29] -= (fwd_rates[180] - rev_rates[180]);
  //sp 7
  sp_rates[6] += (fwd_rates[180] - rev_rates[180]);

  //rxn 181
  //sp 9
  sp_rates[8] -= (fwd_rates[181] - rev_rates[181]);
  //sp 29
  sp_rates[28] += (fwd_rates[181] - rev_rates[181]);
  //sp 5
  sp_rates[4] += (fwd_rates[181] - rev_rates[181]);
  //sp 30
  sp_rates[29] -= (fwd_rates[181] - rev_rates[181]);

  //rxn 182
  //sp 29
  sp_rates[28] += (fwd_rates[182] - rev_rates[182]);
  //sp 5
  sp_rates[4] -= (fwd_rates[182] - rev_rates[182]);
  //sp 30
  sp_rates[29] -= (fwd_rates[182] - rev_rates[182]);
  //sp 8
  sp_rates[7] += (fwd_rates[182] - rev_rates[182]);

  //rxn 183
  //sp 29
  sp_rates[28] += (fwd_rates[183] - rev_rates[183]);
  //sp 12
  sp_rates[11] -= (fwd_rates[183] - rev_rates[183]);
  //sp 13
  sp_rates[12] += (fwd_rates[183] - rev_rates[183]);
  //sp 30
  sp_rates[29] -= (fwd_rates[183] - rev_rates[183]);

  //rxn 184
  //sp 2
  sp_rates[1] -= (fwd_rates[184] - rev_rates[184]);
  //sp 22
  sp_rates[21] += 2.0 * (fwd_rates[184] - rev_rates[184]);
  //sp 31
  sp_rates[30] -= (fwd_rates[184] - rev_rates[184]);

  //rxn 185
  //sp 2
  sp_rates[1] -= (fwd_rates[185] - rev_rates[185]);
  //sp 6
  sp_rates[5] += (fwd_rates[185] - rev_rates[185]);
  //sp 31
  sp_rates[30] -= (fwd_rates[185] - rev_rates[185]);
  //sp 30
  sp_rates[29] += (fwd_rates[185] - rev_rates[185]);

  //rxn 186
  //sp 3
  sp_rates[2] -= (fwd_rates[186] - rev_rates[186]);
  //sp 27
  sp_rates[26] += (fwd_rates[186] - rev_rates[186]);
  //sp 15
  sp_rates[14] += (fwd_rates[186] - rev_rates[186]);
  //sp 31
  sp_rates[30] -= (fwd_rates[186] - rev_rates[186]);

  //rxn 187
  //sp 4
  sp_rates[3] -= (fwd_rates[187] - rev_rates[187]);
  //sp 7
  sp_rates[6] += (fwd_rates[187] - rev_rates[187]);
  //sp 31
  sp_rates[30] -= (fwd_rates[187] - rev_rates[187]);
  //sp 30
  sp_rates[29] += (fwd_rates[187] - rev_rates[187]);

  //rxn 188
  //sp 9
  sp_rates[8] -= (fwd_rates[188] - rev_rates[188]);
  //sp 5
  sp_rates[4] += (fwd_rates[188] - rev_rates[188]);
  //sp 30
  sp_rates[29] += (fwd_rates[188] - rev_rates[188]);
  //sp 31
  sp_rates[30] -= (fwd_rates[188] - rev_rates[188]);

  //rxn 189
  //sp 27
  sp_rates[26] += (fwd_rates[189] - rev_rates[189]);
  //sp 4
  sp_rates[3] += (fwd_rates[189] - rev_rates[189]);
  //sp 5
  sp_rates[4] -= (fwd_rates[189] - rev_rates[189]);
  //sp 15
  sp_rates[14] += (fwd_rates[189] - rev_rates[189]);
  //sp 31
  sp_rates[30] -= (fwd_rates[189] - rev_rates[189]);

  //rxn 190
  //sp 12
  sp_rates[11] -= (fwd_rates[190] - rev_rates[190]);
  //sp 13
  sp_rates[12] += (fwd_rates[190] - rev_rates[190]);
  //sp 30
  sp_rates[29] += (fwd_rates[190] - rev_rates[190]);
  //sp 31
  sp_rates[30] -= (fwd_rates[190] - rev_rates[190]);

  //rxn 191
  //sp 25
  sp_rates[24] += fwd_rates[191];
  //sp 21
  sp_rates[20] += fwd_rates[191];
  //sp 32
  sp_rates[31] = -fwd_rates[191];

  //rxn 192
  //sp 26
  sp_rates[25] += fwd_rates[192];
  //sp 20
  sp_rates[19] += fwd_rates[192];
  //sp 32
  sp_rates[31] -= fwd_rates[192];

  //rxn 193
  //sp 33
  sp_rates[32] = -(fwd_rates[193] - rev_rates[191]) * pres_mod[25];
  //sp 2
  sp_rates[1] -= (fwd_rates[193] - rev_rates[191]) * pres_mod[25];
  //sp 34
  sp_rates[33] = (fwd_rates[193] - rev_rates[191]) * pres_mod[25];

  //rxn 194
  //sp 33
  sp_rates[32] -= (fwd_rates[194] - rev_rates[192]);
  //sp 2
  sp_rates[1] -= (fwd_rates[194] - rev_rates[192]);
  //sp 27
  sp_rates[26] += (fwd_rates[194] - rev_rates[192]);
  //sp 21
  sp_rates[20] += (fwd_rates[194] - rev_rates[192]);

  //rxn 195
  //sp 33
  sp_rates[32] -= (fwd_rates[195] - rev_rates[193]);
  //sp 2
  sp_rates[1] -= (fwd_rates[195] - rev_rates[193]);
  //sp 22
  sp_rates[21] += (fwd_rates[195] - rev_rates[193]);
  //sp 26
  sp_rates[25] += (fwd_rates[195] - rev_rates[193]);

  //rxn 196
  //sp 34
  sp_rates[33] += (fwd_rates[196] - rev_rates[194]);
  //sp 27
  sp_rates[26] -= (fwd_rates[196] - rev_rates[194]);
  //sp 21
  sp_rates[20] -= (fwd_rates[196] - rev_rates[194]);

  //rxn 197
  //sp 2
  sp_rates[1] -= (fwd_rates[197] - rev_rates[195]) * pres_mod[26];
  //sp 35
  sp_rates[34] = -(fwd_rates[197] - rev_rates[195]) * pres_mod[26];
  //sp 36
  sp_rates[35] = (fwd_rates[197] - rev_rates[195]) * pres_mod[26];

  //rxn 198
  //sp 2
  sp_rates[1] -= (fwd_rates[198] - rev_rates[196]);
  //sp 35
  sp_rates[34] -= (fwd_rates[198] - rev_rates[196]);
  //sp 21
  sp_rates[20] += (fwd_rates[198] - rev_rates[196]);
  //sp 31
  sp_rates[30] += (fwd_rates[198] - rev_rates[196]);

  //rxn 199
  //sp 2
  sp_rates[1] -= (fwd_rates[199] - rev_rates[197]);
  //sp 35
  sp_rates[34] -= (fwd_rates[199] - rev_rates[197]);
  //sp 27
  sp_rates[26] += (fwd_rates[199] - rev_rates[197]);
  //sp 26
  sp_rates[25] += (fwd_rates[199] - rev_rates[197]);

  //rxn 200
  //sp 36
  sp_rates[35] += (fwd_rates[200] - rev_rates[198]);
  //sp 21
  sp_rates[20] -= (fwd_rates[200] - rev_rates[198]);
  //sp 31
  sp_rates[30] -= (fwd_rates[200] - rev_rates[198]);

  //rxn 201
  //sp 2
  sp_rates[1] -= (fwd_rates[201] - rev_rates[199]) * pres_mod[27];
  //sp 37
  sp_rates[36] = -(fwd_rates[201] - rev_rates[199]) * pres_mod[27];
  //sp 38
  sp_rates[37] = (fwd_rates[201] - rev_rates[199]) * pres_mod[27];

  //rxn 202
  //sp 2
  sp_rates[1] -= (fwd_rates[202] - rev_rates[200]);
  //sp 37
  sp_rates[36] -= (fwd_rates[202] - rev_rates[200]);
  //sp 21
  sp_rates[20] += (fwd_rates[202] - rev_rates[200]);
  //sp 34
  sp_rates[33] += (fwd_rates[202] - rev_rates[200]);

  //rxn 203
  //sp 2
  sp_rates[1] -= (fwd_rates[203] - rev_rates[201]);
  //sp 37
  sp_rates[36] -= (fwd_rates[203] - rev_rates[201]);
  //sp 31
  sp_rates[30] += (fwd_rates[203] - rev_rates[201]);
  //sp 26
  sp_rates[25] += (fwd_rates[203] - rev_rates[201]);

  //rxn 204
  //sp 34
  sp_rates[33] -= (fwd_rates[204] - rev_rates[202]);
  //sp 21
  sp_rates[20] -= (fwd_rates[204] - rev_rates[202]);
  //sp 38
  sp_rates[37] += (fwd_rates[204] - rev_rates[202]);

  //rxn 205
  //sp 2
  sp_rates[1] -= (fwd_rates[205] - rev_rates[203]) * pres_mod[28];
  //sp 39
  sp_rates[38] = -(fwd_rates[205] - rev_rates[203]) * pres_mod[28];
  //sp 40
  sp_rates[39] = (fwd_rates[205] - rev_rates[203]) * pres_mod[28];

  //rxn 206
  //sp 2
  sp_rates[1] -= (fwd_rates[206] - rev_rates[204]);
  //sp 36
  sp_rates[35] += (fwd_rates[206] - rev_rates[204]);
  //sp 21
  sp_rates[20] += (fwd_rates[206] - rev_rates[204]);
  //sp 39
  sp_rates[38] -= (fwd_rates[206] - rev_rates[204]);

  //rxn 207
  //sp 2
  sp_rates[1] -= (fwd_rates[207] - rev_rates[205]);
  //sp 34
  sp_rates[33] += (fwd_rates[207] - rev_rates[205]);
  //sp 39
  sp_rates[38] -= (fwd_rates[207] - rev_rates[205]);
  //sp 26
  sp_rates[25] += (fwd_rates[207] - rev_rates[205]);

  //rxn 208
  //sp 36
  sp_rates[35] -= (fwd_rates[208] - rev_rates[206]);
  //sp 21
  sp_rates[20] -= (fwd_rates[208] - rev_rates[206]);
  //sp 40
  sp_rates[39] += (fwd_rates[208] - rev_rates[206]);

  //rxn 209
  //sp 41
  sp_rates[40] = -(fwd_rates[209] - rev_rates[207]) * pres_mod[29];
  //sp 2
  sp_rates[1] -= (fwd_rates[209] - rev_rates[207]) * pres_mod[29];
  //sp 42
  sp_rates[41] = (fwd_rates[209] - rev_rates[207]) * pres_mod[29];

  //rxn 210
  //sp 41
  sp_rates[40] -= (fwd_rates[210] - rev_rates[208]);
  //sp 2
  sp_rates[1] -= (fwd_rates[210] - rev_rates[208]);
  //sp 21
  sp_rates[20] += (fwd_rates[210] - rev_rates[208]);
  //sp 38
  sp_rates[37] += (fwd_rates[210] - rev_rates[208]);

  //rxn 211
  //sp 41
  sp_rates[40] -= (fwd_rates[211] - rev_rates[209]);
  //sp 2
  sp_rates[1] -= (fwd_rates[211] - rev_rates[209]);
  //sp 36
  sp_rates[35] += (fwd_rates[211] - rev_rates[209]);
  //sp 26
  sp_rates[25] += (fwd_rates[211] - rev_rates[209]);

  //rxn 212
  //sp 42
  sp_rates[41] += (fwd_rates[212] - rev_rates[210]);
  //sp 21
  sp_rates[20] -= (fwd_rates[212] - rev_rates[210]);
  //sp 38
  sp_rates[37] -= (fwd_rates[212] - rev_rates[210]);

  //rxn 213
  //sp 2
  sp_rates[1] -= (fwd_rates[213] - rev_rates[211]) * pres_mod[30];
  //sp 43
  sp_rates[42] = -(fwd_rates[213] - rev_rates[211]) * pres_mod[30];
  //sp 44
  sp_rates[43] = (fwd_rates[213] - rev_rates[211]) * pres_mod[30];

  //rxn 214
  //sp 2
  sp_rates[1] -= (fwd_rates[214] - rev_rates[212]);
  //sp 43
  sp_rates[42] -= (fwd_rates[214] - rev_rates[212]);
  //sp 21
  sp_rates[20] += (fwd_rates[214] - rev_rates[212]);
  //sp 40
  sp_rates[39] += (fwd_rates[214] - rev_rates[212]);

  //rxn 215
  //sp 2
  sp_rates[1] -= (fwd_rates[215] - rev_rates[213]);
  //sp 43
  sp_rates[42] -= (fwd_rates[215] - rev_rates[213]);
  //sp 38
  sp_rates[37] += (fwd_rates[215] - rev_rates[213]);
  //sp 26
  sp_rates[25] += (fwd_rates[215] - rev_rates[213]);

  //rxn 216
  //sp 44
  sp_rates[43] += (fwd_rates[216] - rev_rates[214]);
  //sp 21
  sp_rates[20] -= (fwd_rates[216] - rev_rates[214]);
  //sp 40
  sp_rates[39] -= (fwd_rates[216] - rev_rates[214]);

  //rxn 217
  //sp 49
  sp_rates[48] = -(fwd_rates[217] - rev_rates[215]);
  //sp 38
  sp_rates[37] += (fwd_rates[217] - rev_rates[215]);
  //sp 32
  sp_rates[31] += (fwd_rates[217] - rev_rates[215]);

  //rxn 218
  //sp 44
  sp_rates[43] -= (fwd_rates[218] - rev_rates[216]);
  //sp 21
  sp_rates[20] -= (fwd_rates[218] - rev_rates[216]);
  //sp 46
  sp_rates[45] = (fwd_rates[218] - rev_rates[216]);

  //rxn 219
  //sp 46
  sp_rates[45] -= (fwd_rates[219] - rev_rates[217]);
  //sp 48
  sp_rates[47] = (fwd_rates[219] - rev_rates[217]);

  //rxn 220
  //sp 26
  sp_rates[25] -= (fwd_rates[220] - rev_rates[218]);
  //sp 47
  sp_rates[46] = (fwd_rates[220] - rev_rates[218]);
  //sp 42
  sp_rates[41] -= (fwd_rates[220] - rev_rates[218]);

  //rxn 221
  //sp 30
  sp_rates[29] -= (fwd_rates[221] - rev_rates[219]);
  //sp 47
  sp_rates[46] += (fwd_rates[221] - rev_rates[219]);
  //sp 40
  sp_rates[39] -= (fwd_rates[221] - rev_rates[219]);

  //rxn 222
  //sp 33
  sp_rates[32] -= (fwd_rates[222] - rev_rates[220]);
  //sp 38
  sp_rates[37] -= (fwd_rates[222] - rev_rates[220]);
  //sp 48
  sp_rates[47] += (fwd_rates[222] - rev_rates[220]);

  //rxn 223
  //sp 43
  sp_rates[42] -= (fwd_rates[223] - rev_rates[221]);
  //sp 22
  sp_rates[21] -= (fwd_rates[223] - rev_rates[221]);
  //sp 48
  sp_rates[47] += (fwd_rates[223] - rev_rates[221]);

  //rxn 224
  //sp 35
  sp_rates[34] -= (fwd_rates[224] - rev_rates[222]);
  //sp 36
  sp_rates[35] -= (fwd_rates[224] - rev_rates[222]);
  //sp 48
  sp_rates[47] += (fwd_rates[224] - rev_rates[222]);

  //rxn 225
  //sp 41
  sp_rates[40] -= (fwd_rates[225] - rev_rates[223]);
  //sp 27
  sp_rates[26] -= (fwd_rates[225] - rev_rates[223]);
  //sp 48
  sp_rates[47] += (fwd_rates[225] - rev_rates[223]);

  //rxn 226
  //sp 34
  sp_rates[33] -= (fwd_rates[226] - rev_rates[224]);
  //sp 37
  sp_rates[36] -= (fwd_rates[226] - rev_rates[224]);
  //sp 48
  sp_rates[47] += (fwd_rates[226] - rev_rates[224]);

  //rxn 227
  //sp 39
  sp_rates[38] -= (fwd_rates[227] - rev_rates[225]);
  //sp 31
  sp_rates[30] -= (fwd_rates[227] - rev_rates[225]);
  //sp 48
  sp_rates[47] += (fwd_rates[227] - rev_rates[225]);

  //rxn 228
  //sp 44
  sp_rates[43] -= (fwd_rates[228] - rev_rates[226]);
  //sp 45
  sp_rates[44] = (fwd_rates[228] - rev_rates[226]);
  //sp 22
  sp_rates[21] -= (fwd_rates[228] - rev_rates[226]);

  //rxn 229
  //sp 42
  sp_rates[41] -= (fwd_rates[229] - rev_rates[227]);
  //sp 27
  sp_rates[26] -= (fwd_rates[229] - rev_rates[227]);
  //sp 45
  sp_rates[44] += (fwd_rates[229] - rev_rates[227]);

  //rxn 230
  //sp 45
  sp_rates[44] += (fwd_rates[230] - rev_rates[228]);
  //sp 31
  sp_rates[30] -= (fwd_rates[230] - rev_rates[228]);
  //sp 40
  sp_rates[39] -= (fwd_rates[230] - rev_rates[228]);

  //rxn 231
  //sp 34
  sp_rates[33] -= (fwd_rates[231] - rev_rates[229]);
  //sp 45
  sp_rates[44] += (fwd_rates[231] - rev_rates[229]);
  //sp 38
  sp_rates[37] -= (fwd_rates[231] - rev_rates[229]);

  //rxn 232
  //sp 36
  sp_rates[35] -= 2.0 * (fwd_rates[232] - rev_rates[230]);
  //sp 45
  sp_rates[44] += (fwd_rates[232] - rev_rates[230]);

  //rxn 233
  //sp 2
  sp_rates[1] -= (fwd_rates[233] - rev_rates[231]);
  //sp 45
  sp_rates[44] -= (fwd_rates[233] - rev_rates[231]);
  //sp 6
  sp_rates[5] += (fwd_rates[233] - rev_rates[231]);
  //sp 46
  sp_rates[45] += (fwd_rates[233] - rev_rates[231]);

  //rxn 234
  //sp 2
  sp_rates[1] -= (fwd_rates[234] - rev_rates[232]);
  //sp 45
  sp_rates[44] -= (fwd_rates[234] - rev_rates[232]);
  //sp 6
  sp_rates[5] += (fwd_rates[234] - rev_rates[232]);
  //sp 47
  sp_rates[46] += (fwd_rates[234] - rev_rates[232]);

  //rxn 235
  //sp 2
  sp_rates[1] -= (fwd_rates[235] - rev_rates[233]);
  //sp 45
  sp_rates[44] -= (fwd_rates[235] - rev_rates[233]);
  //sp 6
  sp_rates[5] += (fwd_rates[235] - rev_rates[233]);
  //sp 48
  sp_rates[47] += (fwd_rates[235] - rev_rates[233]);

  //rxn 236
  //sp 3
  sp_rates[2] -= (fwd_rates[236] - rev_rates[234]);
  //sp 4
  sp_rates[3] += (fwd_rates[236] - rev_rates[234]);
  //sp 45
  sp_rates[44] -= (fwd_rates[236] - rev_rates[234]);
  //sp 46
  sp_rates[45] += (fwd_rates[236] - rev_rates[234]);

  //rxn 237
  //sp 3
  sp_rates[2] -= (fwd_rates[237] - rev_rates[235]);
  //sp 4
  sp_rates[3] += (fwd_rates[237] - rev_rates[235]);
  //sp 45
  sp_rates[44] -= (fwd_rates[237] - rev_rates[235]);
  //sp 47
  sp_rates[46] += (fwd_rates[237] - rev_rates[235]);

  //rxn 238
  //sp 3
  sp_rates[2] -= (fwd_rates[238] - rev_rates[236]);
  //sp 4
  sp_rates[3] += (fwd_rates[238] - rev_rates[236]);
  //sp 45
  sp_rates[44] -= (fwd_rates[238] - rev_rates[236]);
  //sp 48
  sp_rates[47] += (fwd_rates[238] - rev_rates[236]);

  //rxn 239
  //sp 4
  sp_rates[3] -= (fwd_rates[239] - rev_rates[237]);
  //sp 45
  sp_rates[44] -= (fwd_rates[239] - rev_rates[237]);
  //sp 46
  sp_rates[45] += (fwd_rates[239] - rev_rates[237]);
  //sp 7
  sp_rates[6] += (fwd_rates[239] - rev_rates[237]);

  //rxn 240
  //sp 4
  sp_rates[3] -= (fwd_rates[240] - rev_rates[238]);
  //sp 45
  sp_rates[44] -= (fwd_rates[240] - rev_rates[238]);
  //sp 47
  sp_rates[46] += (fwd_rates[240] - rev_rates[238]);
  //sp 7
  sp_rates[6] += (fwd_rates[240] - rev_rates[238]);

  //rxn 241
  //sp 4
  sp_rates[3] -= (fwd_rates[241] - rev_rates[239]);
  //sp 45
  sp_rates[44] -= (fwd_rates[241] - rev_rates[239]);
  //sp 7
  sp_rates[6] += (fwd_rates[241] - rev_rates[239]);
  //sp 48
  sp_rates[47] += (fwd_rates[241] - rev_rates[239]);

  //rxn 242
  //sp 9
  sp_rates[8] -= (fwd_rates[242] - rev_rates[240]);
  //sp 5
  sp_rates[4] += (fwd_rates[242] - rev_rates[240]);
  //sp 45
  sp_rates[44] -= (fwd_rates[242] - rev_rates[240]);
  //sp 46
  sp_rates[45] += (fwd_rates[242] - rev_rates[240]);

  //rxn 243
  //sp 9
  sp_rates[8] -= (fwd_rates[243] - rev_rates[241]);
  //sp 5
  sp_rates[4] += (fwd_rates[243] - rev_rates[241]);
  //sp 45
  sp_rates[44] -= (fwd_rates[243] - rev_rates[241]);
  //sp 47
  sp_rates[46] += (fwd_rates[243] - rev_rates[241]);

  //rxn 244
  //sp 9
  sp_rates[8] -= (fwd_rates[244] - rev_rates[242]);
  //sp 5
  sp_rates[4] += (fwd_rates[244] - rev_rates[242]);
  //sp 45
  sp_rates[44] -= (fwd_rates[244] - rev_rates[242]);
  //sp 48
  sp_rates[47] += (fwd_rates[244] - rev_rates[242]);

  //rxn 245
  //sp 5
  sp_rates[4] -= (fwd_rates[245] - rev_rates[243]);
  //sp 45
  sp_rates[44] -= (fwd_rates[245] - rev_rates[243]);
  //sp 46
  sp_rates[45] += (fwd_rates[245] - rev_rates[243]);
  //sp 8
  sp_rates[7] += (fwd_rates[245] - rev_rates[243]);

  //rxn 246
  //sp 5
  sp_rates[4] -= (fwd_rates[246] - rev_rates[244]);
  //sp 45
  sp_rates[44] -= (fwd_rates[246] - rev_rates[244]);
  //sp 47
  sp_rates[46] += (fwd_rates[246] - rev_rates[244]);
  //sp 8
  sp_rates[7] += (fwd_rates[246] - rev_rates[244]);

  //rxn 247
  //sp 5
  sp_rates[4] -= (fwd_rates[247] - rev_rates[245]);
  //sp 48
  sp_rates[47] += (fwd_rates[247] - rev_rates[245]);
  //sp 45
  sp_rates[44] -= (fwd_rates[247] - rev_rates[245]);
  //sp 8
  sp_rates[7] += (fwd_rates[247] - rev_rates[245]);

  //rxn 248
  //sp 13
  sp_rates[12] += (fwd_rates[248] - rev_rates[246]);
  //sp 12
  sp_rates[11] -= (fwd_rates[248] - rev_rates[246]);
  //sp 45
  sp_rates[44] -= (fwd_rates[248] - rev_rates[246]);
  //sp 46
  sp_rates[45] += (fwd_rates[248] - rev_rates[246]);

  //rxn 249
  //sp 13
  sp_rates[12] += (fwd_rates[249] - rev_rates[247]);
  //sp 12
  sp_rates[11] -= (fwd_rates[249] - rev_rates[247]);
  //sp 45
  sp_rates[44] -= (fwd_rates[249] - rev_rates[247]);
  //sp 47
  sp_rates[46] += (fwd_rates[249] - rev_rates[247]);

  //rxn 250
  //sp 13
  sp_rates[12] += (fwd_rates[250] - rev_rates[248]);
  //sp 12
  sp_rates[11] -= (fwd_rates[250] - rev_rates[248]);
  //sp 45
  sp_rates[44] -= (fwd_rates[250] - rev_rates[248]);
  //sp 48
  sp_rates[47] += (fwd_rates[250] - rev_rates[248]);

  //rxn 251
  //sp 9
  sp_rates[8] -= fwd_rates[251];
  //sp 50
  sp_rates[49] = fwd_rates[251];
  //sp 46
  sp_rates[45] -= fwd_rates[251];

  //rxn 252
  //sp 9
  sp_rates[8] += fwd_rates[252];
  //sp 50
  sp_rates[49] -= fwd_rates[252];
  //sp 46
  sp_rates[45] += fwd_rates[252];

  //rxn 253
  //sp 9
  sp_rates[8] -= fwd_rates[253];
  //sp 50
  sp_rates[49] += fwd_rates[253];
  //sp 47
  sp_rates[46] -= fwd_rates[253];

  //rxn 254
  //sp 9
  sp_rates[8] += fwd_rates[254];
  //sp 50
  sp_rates[49] -= fwd_rates[254];
  //sp 47
  sp_rates[46] += fwd_rates[254];

  //rxn 255
  //sp 9
  sp_rates[8] -= fwd_rates[255];
  //sp 50
  sp_rates[49] += fwd_rates[255];
  //sp 48
  sp_rates[47] -= fwd_rates[255];

  //rxn 256
  //sp 9
  sp_rates[8] += fwd_rates[256];
  //sp 50
  sp_rates[49] -= fwd_rates[256];
  //sp 48
  sp_rates[47] += fwd_rates[256];

  //rxn 257
  //sp 50
  sp_rates[49] -= fwd_rates[257];
  //sp 51
  sp_rates[50] = fwd_rates[257];

  //rxn 258
  //sp 50
  sp_rates[49] += fwd_rates[258];
  //sp 51
  sp_rates[50] -= fwd_rates[258];

  //rxn 259
  //sp 9
  sp_rates[8] -= fwd_rates[259];
  //sp 49
  sp_rates[48] += fwd_rates[259];
  //sp 5
  sp_rates[4] += fwd_rates[259];
  //sp 46
  sp_rates[45] -= fwd_rates[259];

  //rxn 260
  //sp 49
  sp_rates[48] -= fwd_rates[260];
  //sp 9
  sp_rates[8] += fwd_rates[260];
  //sp 5
  sp_rates[4] -= fwd_rates[260];
  //sp 46
  sp_rates[45] += fwd_rates[260];

  //rxn 261
  //sp 9
  sp_rates[8] -= fwd_rates[261];
  //sp 49
  sp_rates[48] += fwd_rates[261];
  //sp 5
  sp_rates[4] += fwd_rates[261];
  //sp 47
  sp_rates[46] -= fwd_rates[261];

  //rxn 262
  //sp 49
  sp_rates[48] -= fwd_rates[262];
  //sp 9
  sp_rates[8] += fwd_rates[262];
  //sp 5
  sp_rates[4] -= fwd_rates[262];
  //sp 47
  sp_rates[46] += fwd_rates[262];

  //rxn 263
  //sp 9
  sp_rates[8] -= fwd_rates[263];
  //sp 49
  sp_rates[48] += fwd_rates[263];
  //sp 5
  sp_rates[4] += fwd_rates[263];
  //sp 48
  sp_rates[47] -= fwd_rates[263];

  //rxn 264
  //sp 49
  sp_rates[48] -= fwd_rates[264];
  //sp 9
  sp_rates[8] += fwd_rates[264];
  //sp 5
  sp_rates[4] -= fwd_rates[264];
  //sp 48
  sp_rates[47] += fwd_rates[264];

  //rxn 265
  //sp 9
  sp_rates[8] -= fwd_rates[265];
  //sp 51
  sp_rates[50] -= fwd_rates[265];
  //sp 52
  sp_rates[51] = fwd_rates[265];

  //rxn 266
  //sp 9
  sp_rates[8] += fwd_rates[266];
  //sp 51
  sp_rates[50] += fwd_rates[266];
  //sp 52
  sp_rates[51] -= fwd_rates[266];

  //rxn 267
  //sp 4
  sp_rates[3] += (fwd_rates[267] - rev_rates[249]);
  //sp 52
  sp_rates[51] -= (fwd_rates[267] - rev_rates[249]);
  //sp 53
  sp_rates[52] = (fwd_rates[267] - rev_rates[249]);

  //rxn 268
  //sp 21
  sp_rates[20] += 3.0 * fwd_rates[268];
  //sp 4
  sp_rates[3] += fwd_rates[268];
  //sp 53
  sp_rates[52] -= fwd_rates[268];
  //sp 22
  sp_rates[21] += fwd_rates[268];
  //sp 24
  sp_rates[23] += 2.0 * fwd_rates[268];

  //sp 1
  sp_rates[0] = 0.0;
  //sp 0
  (*dy_N) = 0.0;
} // end eval_spec_rates

