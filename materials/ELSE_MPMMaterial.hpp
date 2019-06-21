// Base Class of an MPM - Material in ELSE
#ifndef _ELSE_MPM_MATERIAL_HPP_
#define _ELSE_MPM_MATERIAL_HPP_

#include <string>
#include <map>
#include <array>

/*
The Material Base Class
  - The Material holds its properties, as they are constant for all particles this material is called by.
  - The Material provides methods to add a Material Parameter. Their particluar handling and usage is implemented in the child classes.
  - The Material needs a name at its creation.

  Most important interface is the getCauchyStress function:
    return : 0 -> no error
    return : <0 -> error code e.g divergence of local procedure etc...

  Following ELSE standard notation for tensor - vector transition:
    S[6] = [S_11, S_12, S_13, S_22, F_23, S_33]
    F[9] = [F_11, F_12, F_13, F_21, F_22, F_23, F_31, F_32, F_33]
*/

namespace ELSE{
namespace MPM{

class Material {

public:
  // Constructor
  Material(std::string MaterialName);
  // Destructor
  ~Material();

  // Adding Material Parameter
  void addMaterialParameter(std::string Name, int     Parameter);
  void addMaterialParameter(std::string Name, double  Parameter);
  void addMaterialParameter(std::string Name, bool    Parameter);

  // Ask for Value of a Material Parameter
  void getMaterialParameter(std::string KEY, int    &Parameter);
  void getMaterialParameter(std::string KEY, double &Parameter);
  void getMaterialParameter(std::string KEY, bool   &Parameter);

  // getting Cauchy Stressses
  virtual int getCauchyStress(  std::array<double, 6> &CauchyStress,
                                std::array<double, 9> &DeformationGradient,
                                std::map<std::string, double> &MaterialHistory,
                                std::map<std::string, int>    &IntegerMaterialIO,
                                std::map<std::string, double> &DoubleMaterialIO);

private:
  // Material Name
  std::string Name;

  // Materal Parameter container
  std::map<std::string, int     > IntegerParameter;
  std::map<std::string, double  > DoubleParameter;
  std::map<std::string, bool    > BoolParameter;

};

}
}


#endif
