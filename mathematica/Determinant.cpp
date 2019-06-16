/*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           16 Jun 19 13:28:39 *
**************************************************************
User     : Full professional version
Notebook : helper_functions
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 1       Method: Automatic
Subroutine                      : Determinant size: 353
Total size of Mathematica  code : 353 subexpressions
Total size of C code            : 463 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void Determinant(double P[9],double F[9],double T[9])
{
double v[138];
T[0]=F[0]*P[0]+F[1]*P[1]+F[2]*P[2];
T[1]=F[3]*P[0]+F[4]*P[1]+F[5]*P[2];
T[2]=F[6]*P[0]+F[7]*P[1]+F[8]*P[2];
T[3]=F[0]*P[3]+F[1]*P[4]+F[2]*P[5];
T[4]=F[3]*P[3]+F[4]*P[4]+F[5]*P[5];
T[5]=F[6]*P[3]+F[7]*P[4]+F[8]*P[5];
T[6]=F[0]*P[6]+F[1]*P[7]+F[2]*P[8];
T[7]=F[3]*P[6]+F[4]*P[7]+F[5]*P[8];
T[8]=F[6]*P[6]+F[7]*P[7]+F[8]*P[8];
};
