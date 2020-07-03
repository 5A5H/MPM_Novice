/*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           28 Jul 19 11:04:29 *
**************************************************************
User     : Full professional version
Notebook : helper_functions
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 7       Method: Automatic
Subroutine                      : P2S size: 373
Total size of Mathematica  code : 373 subexpressions
Total size of C code            : 503 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void P2S(double Sig[6],double F[9],double P[9])
{
double v[144];
v[19]=F[2]*(-(F[4]*F[6])+F[3]*F[7])+F[1]*(F[5]*F[6]-F[3]*F[8])+F[0]*(-(F[5]*F[7])+F[4]*F[8]);
Sig[0]=(F[0]*P[0]+F[1]*P[1]+F[2]*P[2])/v[19];
Sig[1]=(F[3]*P[0]+F[4]*P[1]+F[5]*P[2])/v[19];
Sig[2]=(F[6]*P[0]+F[7]*P[1]+F[8]*P[2])/v[19];
Sig[3]=(F[3]*P[3]+F[4]*P[4]+F[5]*P[5])/v[19];
Sig[4]=(F[6]*P[3]+F[7]*P[4]+F[8]*P[5])/v[19];
Sig[5]=(F[6]*P[6]+F[7]*P[7]+F[8]*P[8])/v[19];
};
