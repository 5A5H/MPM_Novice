/*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           7 Jun 19 12:10:01  *
**************************************************************
User     : Full professional version
Notebook : helper_functions
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 2       Method: Automatic
Subroutine                      : VonMisesStress size: 120
Total size of Mathematica  code : 120 subexpressions
Total size of C code            : 363 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void VonMisesStress(double Sig[3][3],double (*SigMises))
{
double v[125];
v[11]=(-Sig[0][0]-Sig[1][1]-Sig[2][2])/3e0;
(*SigMises)=sqrt(0.15e1*(2e0*Sig[0][1]*Sig[1][0]+2e0*Sig[0][2]*Sig[2][0]+2e0*Sig[1][2]*Sig[2][1]
 +Power(Sig[0][0]+v[11],2)+Power(Sig[1][1]+v[11],2)+Power(Sig[2][2]+v[11],2)));
};
