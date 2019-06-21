#ifndef _LINEAR_ELASTICITY_A_HPP_
#define _LINEAR_ELASTICITY_A_HPP_

// inlcude the base class
#include <ELSE_MPMMaterial.hpp>

#include <string>
#include <map>
#include <array>
#include <iostream>

/*
Imolementation of a Linear Elasticity Law
Files:
      materials/LinearElasticity_A.hpp
      materials/LinearElasticity_A.cpp

  Following ELSE standard notation for tensor - vector transition:
    S[6] = [S_11, S_12, S_13, S_22, F_23, S_33]
    F[9] = [F_11, F_12, F_13, F_21, F_22, F_23, F_31, F_32, F_33]
*/

namespace ELSE{
namespace MPM{

class LinearElasticity_A : public Material {

public:
  // Constructor
  LinearElasticity_A(std::string MaterialName);
  // Destructor
  ~LinearElasticity_A();

  // getting Cauchy Stressses
  int getCauchyStress(  std::array<double, 6> &CauchyStress,
                        std::array<double, 9> &DeformationGradient,
                        std::map<std::string, double> &MaterialHistory,
                        std::map<std::string, int>    &IntegerMaterialIO,
                        std::map<std::string, double> &DoubleMaterialIO);

private:
  // Check wether all needed Material Parameters are present.
  int checkMaterialParameters();

};

}
}


#endif
