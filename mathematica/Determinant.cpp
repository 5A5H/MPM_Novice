/*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           16 Jun 19 22:57:20 *
**************************************************************
User     : Full professional version
Notebook : helper_functions
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 1       Method: Automatic
Subroutine                      : Determinant size: 113
Total size of Mathematica  code : 113 subexpressions
Total size of C code            : 279 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void Determinant(double F[3][3],double (*detF))
{
double v[121];
(*detF)=F[0][2]*(-(F[1][1]*F[2][0])+F[1][0]*F[2][1])-F[0][1]*(-(F[1][2]*F[2][0])+F[1][0]*F[2][2])
 +F[0][0]*(-(F[1][2]*F[2][1])+F[1][1]*F[2][2]);
};
