/*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           26 May 19 17:13:05 *
**************************************************************
User     : Full professional version
Notebook : mate_planestrain_linearelasticity_ctest
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 9       Method: Automatic
Subroutine                      : SmallStrainHookePlaneStress2D size: 233
Total size of Mathematica  code : 233 subexpressions
Total size of C code            : 651 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void SmallStrainHookePlaneStress2D(double MaterialData[2],double Eps[3]
     ,double Exp[16])
{
double v[144];
v[34]=Eps[0]+Eps[1];
v[32]=MaterialData[0]/(1e0+MaterialData[1]);
v[3]=(MaterialData[1]*v[32])/(1e0-2e0*MaterialData[1]);
v[33]=-(v[3]/(v[3]+v[32]));
v[28]=v[3]*(1e0+v[33]);
v[30]=v[28]+v[32];
v[20]=v[3]*(v[34]+v[33]*v[34]);
Exp[0]=v[20]+Eps[0]*v[32];
Exp[1]=v[20]+Eps[1]*v[32];
Exp[2]=0e0;
Exp[3]=Eps[2]*v[32];
Exp[4]=v[30];
Exp[5]=v[28];
Exp[6]=0e0;
Exp[7]=v[28];
Exp[8]=v[30];
Exp[9]=0e0;
Exp[10]=0e0;
Exp[11]=0e0;
Exp[12]=0e0;
Exp[13]=0e0;
Exp[14]=0e0;
Exp[15]=v[32];
};
