/*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           15 Jun 19 22:57:34 *
**************************************************************
User     : Full professional version
Notebook : mate_planestrain_linearelasticity_ctest
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 12      Method: Automatic
Subroutine                      : SmallStrainHookePlaneStress2D size: 243
Total size of Mathematica  code : 243 subexpressions
Total size of C code            : 756 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void SmallStrainHookePlaneStress2D(double MaterialData[2],double Eps[3]
     ,double Exp[16])
{
double v[149];
v[39]=Eps[0]+Eps[1];
bool debugaceroutine=true;
v[37]=MaterialData[0]/(1e0+MaterialData[1]);
v[4]=(MaterialData[1]*v[37])/(1e0-2e0*MaterialData[1]);
v[38]=-(v[4]/(v[37]+v[4]));
v[31]=(1e0+v[38])*v[4];
v[33]=v[31]+v[37];
v[21]=(v[39]+v[38]*v[39])*v[4];
v[18]=v[21]+Eps[0]*v[37];
Exp[0]=v[18];
Exp[1]=v[21]+Eps[1]*v[37];
Exp[2]=0e0;
Exp[3]=Eps[2]*v[37];
if (debugaceroutine) std::cout << "sigt: " <<v[18] << std::endl;
Exp[4]=v[33];
Exp[5]=v[31];
Exp[6]=0e0;
Exp[7]=v[31];
Exp[8]=v[33];
Exp[9]=0e0;
Exp[10]=0e0;
Exp[11]=0e0;
Exp[12]=0e0;
Exp[13]=0e0;
Exp[14]=0e0;
Exp[15]=v[37];
};
