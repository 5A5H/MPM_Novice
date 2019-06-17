/*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           16 Jun 19 23:15:44 *
**************************************************************
User     : Full professional version
Notebook : helper_functions
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 1       Method: Automatic
Subroutine                      : Tensor_ijkj size: 353
Total size of Mathematica  code : 353 subexpressions
Total size of C code            : 463 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void Tensor_ijkj(double A[9],double B[9],double C[9])
{
double v[138];
C[0]=A[0]*B[0]+A[1]*B[1]+A[2]*B[2];
C[1]=A[0]*B[3]+A[1]*B[4]+A[2]*B[5];
C[2]=A[0]*B[6]+A[1]*B[7]+A[2]*B[8];
C[3]=A[3]*B[0]+A[4]*B[1]+A[5]*B[2];
C[4]=A[3]*B[3]+A[4]*B[4]+A[5]*B[5];
C[5]=A[3]*B[6]+A[4]*B[7]+A[5]*B[8];
C[6]=A[6]*B[0]+A[7]*B[1]+A[8]*B[2];
C[7]=A[6]*B[3]+A[7]*B[4]+A[8]*B[5];
C[8]=A[6]*B[6]+A[7]*B[7]+A[8]*B[8];
};
