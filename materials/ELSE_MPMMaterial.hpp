// Base Class of an MPM - Material in ELSE
#ifndef _ELSE_MPM_MATERIAL_HPP_
#define _ELSE_MPM_MATERIAL_HPP_

/*
The Material Base Class
*/

namespace ELSE{
namespace MPM{

class Material {

public:
  // Constructor
  Material();
  // Destructor
  ~Material();

  // Adding Material Parameter
  void addMaterialParameter(std::string Name, int     Parameter);
  void addMaterialParameter(std::string Name, double  Parameter);
  void addMaterialParameter(std::string Name, bool    Parameter);

private:
};

// Constructor implementation
Material::Material(){}

// Destructor implementation
Material::~Material(){}


}
}


#endif
