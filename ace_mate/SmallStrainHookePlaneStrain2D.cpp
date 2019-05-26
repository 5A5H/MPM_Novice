/*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           26 May 19 17:12:41 *
**************************************************************
User     : Full professional version
Notebook : mate_planestrain_linearelasticity_ctest
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 6       Method: Automatic
Subroutine                      : SmallStrainHookePlaneStrain2D size: 205
Total size of Mathematica  code : 205 subexpressions
Total size of C code            : 582 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void SmallStrainHookePlaneStrain2D(double MaterialData[2],double Eps[3]
     ,double Exp[16])
{
double v[138];
v[28]=MaterialData[0]/(1e0+MaterialData[1]);
v[12]=(MaterialData[1]*v[28])/(1e0-2e0*MaterialData[1]);
v[20]=(Eps[0]+Eps[1])*v[12];
v[25]=v[12]+v[28];
Exp[0]=v[20]+Eps[0]*v[28];
Exp[1]=v[20]+Eps[1]*v[28];
Exp[2]=v[20];
Exp[3]=Eps[2]*v[28];
Exp[4]=v[25];
Exp[5]=v[12];
Exp[6]=0e0;
Exp[7]=v[12];
Exp[8]=v[25];
Exp[9]=0e0;
Exp[10]=v[12];
Exp[11]=v[12];
Exp[12]=0e0;
Exp[13]=0e0;
Exp[14]=0e0;
Exp[15]=v[28];
};
