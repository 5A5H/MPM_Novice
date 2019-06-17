#ifndef _ELSE_CONTINUUM_MECHANICS_LIBRARY_HPP_
#define _ELSE_CONTINUUM_MECHANICS_LIBRARY_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <math.h>

namespace ELSE {
namespace Conti {

// VonMises stress computation
inline void VonMisesStress(double Sig[3][3],double &SigMises){SigMises=sqrt(0.15e1*(2e0*Sig[0][1]*Sig[1][0]+2e0*Sig[0][2]*Sig[2][0]+2e0*Sig[1][2]*Sig[2][1]+pow(Sig[0][0]+((-Sig[0][0]-Sig[1][1]-Sig[2][2])/3e0),2)+pow(Sig[1][1]+((-Sig[0][0]-Sig[1][1]-Sig[2][2])/3e0),2)+pow(Sig[2][2]+((-Sig[0][0]-Sig[1][1]-Sig[2][2])/3e0),2)));};
inline void VonMisesStress(double Sig[9],double &SigMises){SigMises=sqrt(0.15e1*(2e0*Sig[1]*Sig[3]+2e0*Sig[2]*Sig[6]+2e0*Sig[5]*Sig[7]+pow(Sig[0]+((-Sig[0]-Sig[4]-Sig[8])/3e0),2)+pow(Sig[4]+((-Sig[0]-Sig[4]-Sig[8])/3e0),2)+pow(Sig[8]+((-Sig[0]-Sig[4]-Sig[8])/3e0),2)));};

// Determinant computation
inline void Det(double F[9],double &detF){detF=F[2]*(-(F[4]*F[6])+F[3]*F[7])-F[1]*(-(F[5]*F[6])+F[3]*F[8])+F[0]*(-(F[5]*F[7])+F[4]*F[8]);};
inline void Det(double F[3][3],double &detF){detF=F[0][2]*(-(F[1][1]*F[2][0])+F[1][0]*F[2][1])-F[0][1]*(-(F[1][2]*F[2][0])+F[1][0]*F[2][2])+F[0][0]*(-(F[1][2]*F[2][1])+F[1][1]*F[2][2]);};

// Tensor operations
inline void TensorProduct(double A[9],double B[9],double C[9],std::string KEY){
  // This function supports tensor products for the "tensor as vector" notation A[9] = [A11, A12, A13, A21, A22, A23, A31, A32, A33]
  // KEY  | Math       | Index
  //------------------------------------
  // ijij | C = A:B    | C = A_ij B_ij
  // ijjk | C = A.B    | C = A_ij B_jk
  // ijkj | C = A.B^T  | C = A_ij B_kj
  // jijk | C = A^T.B  | C = A_ji B_jk
  if (KEY=="ijij") {C[0]=A[0]*B[0]+A[1]*B[1]+A[2]*B[2]+A[3]*B[3]+A[4]*B[4]+A[5]*B[5]+A[6]*B[6]+A[7]*B[7]+A[8]*B[8];return;};
  if (KEY=="ijjk") {C[0]=A[0]*B[0]+A[1]*B[3]+A[2]*B[6];C[1]=A[0]*B[1]+A[1]*B[4]+A[2]*B[7];C[2]=A[0]*B[2]+A[1]*B[5]+A[2]*B[8];C[3]=A[3]*B[0]+A[4]*B[3]+A[5]*B[6];C[4]=A[3]*B[1]+A[4]*B[4]+A[5]*B[7];C[5]=A[3]*B[2]+A[4]*B[5]+A[5]*B[8];C[6]=A[6]*B[0]+A[7]*B[3]+A[8]*B[6];C[7]=A[6]*B[1]+A[7]*B[4]+A[8]*B[7];C[8]=A[6]*B[2]+A[7]*B[5]+A[8]*B[8];return;};
  if (KEY=="ijkj") {C[0]=A[0]*B[0]+A[1]*B[1]+A[2]*B[2];C[1]=A[0]*B[3]+A[1]*B[4]+A[2]*B[5];C[2]=A[0]*B[6]+A[1]*B[7]+A[2]*B[8];C[3]=A[3]*B[0]+A[4]*B[1]+A[5]*B[2];C[4]=A[3]*B[3]+A[4]*B[4]+A[5]*B[5];C[5]=A[3]*B[6]+A[4]*B[7]+A[5]*B[8];C[6]=A[6]*B[0]+A[7]*B[1]+A[8]*B[2];C[7]=A[6]*B[3]+A[7]*B[4]+A[8]*B[5];C[8]=A[6]*B[6]+A[7]*B[7]+A[8]*B[8];return;};
  if (KEY=="jijk") {C[0]=A[0]*B[0]+A[3]*B[3]+A[6]*B[6];C[1]=A[0]*B[1]+A[3]*B[4]+A[6]*B[7];C[2]=A[0]*B[2]+A[3]*B[5]+A[6]*B[8];C[3]=A[1]*B[0]+A[4]*B[3]+A[7]*B[6];C[4]=A[1]*B[1]+A[4]*B[4]+A[7]*B[7];C[5]=A[1]*B[2]+A[4]*B[5]+A[7]*B[8];C[6]=A[2]*B[0]+A[5]*B[3]+A[8]*B[6];C[7]=A[2]*B[1]+A[5]*B[4]+A[8]*B[7];C[8]=A[2]*B[2]+A[5]*B[5]+A[8]*B[8];return;};
};

} // namespace Conti
} // namspace ELSE

#endif
