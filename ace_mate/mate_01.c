/*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           24 May 19 09:27:47 *
**************************************************************
User     : Full professional version
Notebook : mate_planestrain_linearelasticity_ctest
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 6       Method: Automatic
Subroutine                      : mate01 size: 205
Total size of Mathematica  code : 205 subexpressions
Total size of C code            : 558 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void mate01(double v[138],double MaterialData[2],double Eps[3]
     ,double Exp[16])
{
v[28]=MaterialData[0]/(1e0+MaterialData[1]);
v[12]=(MaterialData[1]*v[28])/(1e0-2e0*MaterialData[1]);
v[20]=(Eps[0]+Eps[2])*v[12];
v[25]=v[12]+v[28];
Exp[0]=v[20]+Eps[0]*v[28];
Exp[1]=Eps[1]*v[28];
Exp[2]=v[20]+Eps[2]*v[28];
Exp[3]=v[20];
Exp[4]=v[25];
Exp[5]=0e0;
Exp[6]=v[12];
Exp[7]=0e0;
Exp[8]=v[28];
Exp[9]=0e0;
Exp[10]=v[12];
Exp[11]=0e0;
Exp[12]=v[25];
Exp[13]=v[12];
Exp[14]=0e0;
Exp[15]=v[12];
};
