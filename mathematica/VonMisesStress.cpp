/*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           16 Jun 19 22:45:43 *
**************************************************************
User     : Full professional version
Notebook : helper_functions
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 2       Method: Automatic
Subroutine                      : VonMisesStress size: 108
Total size of Mathematica  code : 108 subexpressions
Total size of C code            : 324 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void VonMisesStress(double Sig[9],double (*SigMises))
{
double v[125];
v[11]=(-Sig[0]-Sig[4]-Sig[8])/3e0;
(*SigMises)=sqrt(0.15e1*(2e0*Sig[1]*Sig[3]+2e0*Sig[2]*Sig[6]+2e0*Sig[5]*Sig[7]+Power(Sig[0]+v[11],2
 )+Power(Sig[4]+v[11],2)+Power(Sig[8]+v[11],2)));
};
