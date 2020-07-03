/*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           23 Jul 19 10:39:59 *
**************************************************************
User     : Full professional version
Notebook : helper_functions
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 4       Method: Automatic
Subroutine                      : updateF size: 425
Total size of Mathematica  code : 425 subexpressions
Total size of C code            : 647 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void updateF(double Fn[9],double L[9],double (*dt),double F[9])
{
double v[133];
v[20]=1e0+(*dt)*L[0];
v[21]=1e0+(*dt)*L[4];
v[22]=1e0+(*dt)*L[8];
F[0]=(*dt)*(Fn[3]*L[1]+Fn[6]*L[2])+Fn[0]*v[20];
F[1]=(*dt)*(Fn[4]*L[1]+Fn[7]*L[2])+Fn[1]*v[20];
F[2]=(*dt)*(Fn[5]*L[1]+Fn[8]*L[2])+Fn[2]*v[20];
F[3]=(*dt)*(Fn[0]*L[3]+Fn[6]*L[5])+Fn[3]*v[21];
F[4]=(*dt)*(Fn[1]*L[3]+Fn[7]*L[5])+Fn[4]*v[21];
F[5]=(*dt)*(Fn[2]*L[3]+Fn[8]*L[5])+Fn[5]*v[21];
F[6]=(*dt)*(Fn[0]*L[6]+Fn[3]*L[7])+Fn[6]*v[22];
F[7]=(*dt)*(Fn[1]*L[6]+Fn[4]*L[7])+Fn[7]*v[22];
F[8]=(*dt)*(Fn[2]*L[6]+Fn[5]*L[7])+Fn[8]*v[22];
};
