// Definition of the MPM Material Factory
#ifndef _ELSE_MPM_MATERIAL_FACTORY_HPP_
#define _ELSE_MPM_MATERIAL_FACTORY_HPP_

#include <string>
#include <map>
#include <array>

// Inlcude Material Base Class
#include <ELSE_MPMMaterial.hpp>

// Below add inlcudes for the derived classes
#include <LinearElasticity_A.hpp>

/*
The Material Factory
  - The Factory is called to create material objects.
  - It returns a pointer to these objects.
  - The different material classes need to be implemented as follows:
    * Declaration an Implementation:
                                    materials/Material_A.cpp
                                    materials/Material_A.hpp
    *! Unfortunately by now one has to also add the .cpp to the CMakeList.txt
    * They MUST inherit the MPM Material Base class:
                                    class Material_A: public Material {...}
    * They need to be known within this Factory:
                                    #inlcude <Material_A.hpp>
                      ---to be continued---

  Following ELSE standard notation for tensor - vector transition:
    S[6] = [S_11, S_12, S_13, S_22, F_23, S_33]
    F[9] = [F_11, F_12, F_13, F_21, F_22, F_23, F_31, F_32, F_33]
*/

namespace ELSE{
namespace MPM{

  Material* CreateMaterial(const std::string MaterialKind, const std::string MaterialName){
    // MaterialKey defines the actual type of material object e.g. Material_A
    // MaterialName defines the name that is used for the specific problem e.g. Steel,Alu ..

    if (MaterialKey=="LinearElasticity_A") return new LinearElasticity_A(MaterialName);

    // if function evaluates to here there is no implementation for the requested material key
    std::cerr << "Error: No Implementation for a Material: " << MaterialKind << " found." << std::endl;
    return nullptr;
  }



}
}


#endif
