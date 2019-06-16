#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <math.h>

namespace ELSE {
namespace Geometric {

  inline bool PointInQ4(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
    bool DetailedOutput = false;
    if (DetailedOutput){
    std::cout << "Test X1: " << X1[0] << " ," << X1[1] << " ," << X1[2] << std::endl;
    std::cout << "Test X2: " << X2[0] << " ," << X2[1] << " ," << X2[2] << std::endl;
    std::cout << "Test X3: " << X3[0] << " ," << X3[1] << " ," << X3[2] << std::endl;
    std::cout << "Test X4: " << X4[0] << " ," << X4[1] << " ," << X4[2] << std::endl;
    std::cout << "Test XP: " << XP[0] << " ," << XP[1] << " ," << XP[2] << std::endl;
    }

    // computation of tangents (in element order)
    double t[4][3];
    for(int j=0;j<3;j++){
        t[0][j] = X2[j]-X1[j];
        t[1][j] = X3[j]-X2[j];
        t[2][j] = X4[j]-X3[j];
        t[3][j] = X1[j]-X4[j];
    }
    if (DetailedOutput){
    std::cout << "TangentialVector X1X2: " << t[0][0] << ", " << t[0][1] << ", " << t[0][2] << ", " << std::endl;
    std::cout << "TangentialVector X2X3: " << t[1][0] << ", " << t[1][1] << ", " << t[1][2] << ", " << std::endl;
    std::cout << "TangentialVector X3X4: " << t[2][0] << ", " << t[2][1] << ", " << t[2][2] << ", " << std::endl;
    std::cout << "TangentialVector X4X1: " << t[3][0] << ", " << t[3][1] << ", " << t[3][2] << ", " << std::endl;
    }

    // compute normals
    double n[4][3];
    for(int j=0;j<4;j++){
    n[j][0] = -t[j][1]; n[j][1] = t[j][0]; n[j][2] = 0.0;
    }
    if (DetailedOutput){
    std::cout << "NormalVector X1X2: " << n[0][0] << ", " << n[0][1] << ", " << n[0][2] << ", " << std::endl;
    std::cout << "NormalVector X2X3: " << n[1][0] << ", " << n[1][1] << ", " << n[1][2] << ", " << std::endl;
    std::cout << "NormalVector X3X4: " << n[2][0] << ", " << n[2][1] << ", " << n[2][2] << ", " << std::endl;
    std::cout << "NormalVector X4X1: " << n[3][0] << ", " << n[3][1] << ", " << n[3][2] << ", " << std::endl;
    }
    // compute testvectors
    double tv[4][3];
    for(int j=0;j<3;j++){
        tv[0][j] = XP[j]-X1[j];
        tv[1][j] = XP[j]-X2[j];
        tv[2][j] = XP[j]-X3[j];
        tv[3][j] = XP[j]-X4[j];
    }
    if (DetailedOutput){
    std::cout << "TestVector X1XP: " << tv[0][0] << ", " << tv[0][1] << ", " << tv[0][2] << ", " << std::endl;
    std::cout << "TestVector X2XP: " << tv[1][0] << ", " << tv[1][1] << ", " << tv[1][2] << ", " << std::endl;
    std::cout << "TestVector X3XP: " << tv[2][0] << ", " << tv[2][1] << ", " << tv[2][2] << ", " << std::endl;
    std::cout << "TestVector X4XP: " << tv[3][0] << ", " << tv[3][1] << ", " << tv[3][2] << ", " << std::endl;
    }
    // compute sum of scalar products to comare with abs of scalar products
    double sumabs = 0.0;
    double sum    = 0.0;
    for(int i=0;i<4;i++){
    //std::cout << "ScalarProduct: " << ( tv[i][0]*n[i][0] + tv[i][1]*n[i][1] + tv[i][2]*n[i][2] ) << std::endl;
          sumabs  += abs( tv[i][0]*n[i][0] + tv[i][1]*n[i][1] + tv[i][2]*n[i][2] );
          sum     += ( tv[i][0]*n[i][0] + tv[i][1]*n[i][1] + tv[i][2]*n[i][2] );
    }
    bool XPInside = (sumabs-sum < 10e-10)? true : false;
    if (DetailedOutput){
    std::cout << "sum      : " << sum << std::endl;
    std::cout << "sumabs   : " << sumabs << std::endl;
    if(XPInside) {
        std::cout << "XPInside : True" << std::endl;
      } else {
        std::cout << "XPInside : False" << std::endl;
      }
    }
    return XPInside;
  };


} // namespace Geometric
} // namspace ELSE
