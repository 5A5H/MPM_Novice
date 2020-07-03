/*************************************************************
* AceGen    6.923 MacOSX (19 Apr 19)                         *
*           Co. J. Korelc  2013           27 Jul 19 16:13:53 *
**************************************************************
User     : Full professional version
Notebook : mate_planestrain_linearelasticity_ctest
Evaluation time                 : 1 s     Mode  : Optimal
Number of formulae              : 31      Method: Automatic
Subroutine                      : NeoHookeTypeA size: 1014
Total size of Mathematica  code : 1014 subexpressions
Total size of C code            : 1799 bytes */
#include "sms.h"

/******************* S U B R O U T I N E *********************/
void NeoHookeTypeA(double MaterialData[2],double F[9],double Sig[9]
     ,double Cmat[9][9])
{
double v[173];
v[62]=MaterialData[0]/(1e0+MaterialData[1]);
v[3]=(MaterialData[1]*v[62])/(1e0-2e0*MaterialData[1]);
v[39]=v[62]/4e0;
v[24]=F[2]*(-(F[4]*F[6])+F[3]*F[7])+F[1]*(F[5]*F[6]-F[3]*F[8])+F[0]*(-(F[5]*F[7])+F[4]*F[8]);
v[25]=Power(F[0],2)+Power(F[3],2)+Power(F[6],2);
v[26]=F[0]*F[1]+F[3]*F[4]+F[6]*F[7];
v[35]=(v[26]*v[26]);
v[27]=F[0]*F[2]+F[3]*F[5]+F[6]*F[8];
v[40]=(v[27]*v[27]);
v[37]=2e0*v[26]*v[27];
v[28]=Power(F[1],2)+Power(F[4],2)+Power(F[7],2);
v[29]=F[1]*F[2]+F[4]*F[5]+F[7]*F[8];
v[30]=Power(F[2],2)+Power(F[5],2)+Power(F[8],2);
v[63]=(v[29]*v[29])-v[28]*v[30];
v[34]=(v[3]+(v[3]+v[62])/(v[30]*v[35]-v[29]*v[37]+v[28]*v[40]+v[25]*v[63]))/4e0;
v[36]=v[34]*(v[25]*v[28]-v[35])+v[39];
v[38]=v[34]*(-2e0*v[25]*v[29]+v[37]);
v[41]=v[39]+v[34]*(v[25]*v[30]-v[40]);
v[42]=(-2e0*v[27]*v[28]+2e0*v[26]*v[29])*v[34];
v[43]=2e0*(v[27]*v[29]-v[26]*v[30])*v[34];
v[45]=v[39]-v[34]*v[63];
v[46]=F[2]*v[42]+F[1]*v[43]+2e0*F[0]*v[45];
v[47]=F[2]*v[38]+2e0*F[1]*v[41]+F[0]*v[43];
v[48]=2e0*F[2]*v[36]+F[1]*v[38]+F[0]*v[42];
v[49]=F[5]*v[42]+F[4]*v[43]+2e0*F[3]*v[45];
v[50]=F[5]*v[38]+2e0*F[4]*v[41]+F[3]*v[43];
v[51]=2e0*F[5]*v[36]+F[4]*v[38]+F[3]*v[42];
v[56]=(F[3]*v[46]+F[4]*v[47]+F[5]*v[48])/v[24];
v[57]=(F[6]*v[46]+F[7]*v[47]+F[8]*v[48])/v[24];
v[59]=(F[6]*v[49]+F[7]*v[50]+F[8]*v[51])/v[24];
Sig[0]=(F[0]*v[46]+F[1]*v[47]+F[2]*v[48])/v[24];
Sig[1]=v[56];
Sig[2]=v[57];
Sig[3]=v[56];
Sig[4]=(F[3]*v[49]+F[4]*v[50]+F[5]*v[51])/v[24];
Sig[5]=v[59];
Sig[6]=v[57];
Sig[7]=v[59];
Sig[8]=(F[8]*(2e0*F[8]*v[36]+F[7]*v[38]+F[6]*v[42])+F[7]*(F[8]*v[38]+2e0*F[7]*v[41]+F[6]*v[43])
 +F[6]*(F[8]*v[42]+F[7]*v[43]+2e0*F[6]*v[45]))/v[24];
};
