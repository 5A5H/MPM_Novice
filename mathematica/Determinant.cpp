/*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           16 Jun 19 12:37:11 *
**************************************************************
User     : Full professional version
Notebook : helper_functions
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 1       Method: Automatic
Subroutine                      : Determinant size: 98
Total size of Mathematica  code : 98 subexpressions
Total size of C code            : 229 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void Determinant(double F[9],double (*detF))
{
double v[121];
(*detF)=F[2]*(-(F[4]*F[6])+F[3]*F[7])-F[1]*(-(F[5]*F[6])+F[3]*F[8])+F[0]*(-(F[5]*F[7])+F[4]*F[8]);
};
