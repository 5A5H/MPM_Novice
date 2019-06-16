/*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           15 Jun 19 21:57:02 *
**************************************************************
User     : Full professional version
Notebook : mate_planestrain_linearelasticity_ctest
Evaluation time                 : 1 s     Mode  : Optimal
Number of formulae              : 23      Method: Automatic
Subroutine                      : STVKPlaneStress2D size: 847
Total size of Mathematica  code : 847 subexpressions
Total size of C code            : 1446 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void STVKPlaneStress2D(double MaterialData[2],double F[9],double Sig[9]
     ,double Cmat[9][9])
{
double v[153];
v[42]=MaterialData[0]/(1e0+MaterialData[1]);
v[3]=(MaterialData[1]*v[42])/(1e0-2e0*MaterialData[1]);
v[5]=v[42]/2e0;
v[43]=2e0*(v[3]+v[5]);
v[25]=(2e0*v[5])/(v[3]+v[42]);
v[15]=F[2]*(-(F[4]*F[6])+F[3]*F[7])+F[1]*(F[5]*F[6]-F[3]*F[8])+F[0]*(-(F[5]*F[7])+F[4]*F[8]);
v[16]=(-1e0+Power(F[0],2)+Power(F[3],2)+Power(F[6],2))/2e0;
v[19]=(-1e0+Power(F[1],2)+Power(F[4],2)+Power(F[7],2))/2e0;
v[22]=v[25]*(v[19]*v[3]+v[16]*v[43]);
v[23]=(F[0]*F[1]+F[3]*F[4]+F[6]*F[7])*v[5];
v[24]=(F[0]*F[2]+F[3]*F[5]+F[6]*F[8])*v[5];
v[36]=F[3]*v[22]+F[4]*v[23]+F[5]*v[24];
v[30]=F[0]*v[22]+F[1]*v[23]+F[2]*v[24];
v[27]=v[25]*(v[16]*v[3]+v[19]*v[43]);
v[28]=(F[1]*F[2]+F[4]*F[5]+F[7]*F[8])*v[5];
v[38]=F[3]*v[23]+F[4]*v[27]+F[5]*v[28];
v[37]=F[3]*v[24]+F[4]*v[28];
v[32]=F[0]*v[23]+F[1]*v[27]+F[2]*v[28];
v[31]=F[0]*v[24]+F[1]*v[28];
v[33]=(F[3]*v[30]+F[5]*v[31]+F[4]*v[32])/v[15];
v[34]=(F[6]*v[30]+F[8]*v[31]+F[7]*v[32])/v[15];
v[39]=(F[6]*v[36]+F[8]*v[37]+F[7]*v[38])/v[15];
Sig[0]=(F[0]*v[30]+F[2]*v[31]+F[1]*v[32])/v[15];
Sig[1]=v[33];
Sig[2]=v[34];
Sig[3]=v[33];
Sig[4]=(F[3]*v[36]+F[5]*v[37]+F[4]*v[38])/v[15];
Sig[5]=v[39];
Sig[6]=v[34];
Sig[7]=v[39];
Sig[8]=(F[6]*(F[6]*v[22]+F[7]*v[23]+F[8]*v[24])+F[8]*(F[6]*v[24]+F[7]*v[28])+F[7]*(F[6]*v[23]
 +F[7]*v[27]+F[8]*v[28]))/v[15];
};
