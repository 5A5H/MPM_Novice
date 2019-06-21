#include <ELSE_System.hpp>
#include <ELSE_MPMMaterial.hpp>

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

// Constructor implementation
Material::Material(std::string MaterialName){
  Name = MaterialName;
  if (LogFile.is_open()) LogFile << "Created:\t MPM-Material Object: " << Name << std::endl;
}

// Destructor implementation
Material::~Material(){}

// Dump all available material parameter into offstream
void Material::dumpMaterialParameter(std::ofstream &Outstream){
  Outstream << "Material Parameter of MPM Material: " << Name << std::endl;

  Outstream << "- IntegerParameters\t: " << std::endl;
  for (auto &ParameterEntry : IntegerParameter) {
    Outstream << "\t\t\t \"" << ParameterEntry.first << "\"";
    Outstream << "  ";
    Outstream << ParameterEntry.second;
    Outstream << std::endl;
  }

  Outstream << "- DoubleParameters\t: " << std::endl;
  for (auto &ParameterEntry : DoubleParameter) {
    Outstream << "\t\t\t \"" << ParameterEntry.first << "\"";
    Outstream << "  ";
    Outstream << ParameterEntry.second;
    Outstream << std::endl;
  }

  Outstream << "- BoolParameters\t\t: " << std::endl;
  for (auto &ParameterEntry : BoolParameter) {
    Outstream << "\t\t\t \"" << ParameterEntry.first << "\"";
    Outstream << "  ";
    Outstream << ParameterEntry.second;
    Outstream << std::endl;
  }
};

// Adding Material Parameter (Overloading)
void Material::addMaterialParameter(std::string Name, int     Parameter){IntegerParameter[Name] = Parameter;};
void Material::addMaterialParameter(std::string Name, double  Parameter){DoubleParameter[Name]  = Parameter;};
void Material::addMaterialParameter(std::string Name, bool    Parameter){BoolParameter[Name]    = Parameter;};

// Reading Material Parameter (Overloading)
void Material::getMaterialParameter(std::string KEY, int    &Parameter){Parameter = IntegerParameter.find(KEY)->second; };
void Material::getMaterialParameter(std::string KEY, double &Parameter){Parameter = DoubleParameter.find(KEY)->second;  };
void Material::getMaterialParameter(std::string KEY, bool   &Parameter){Parameter = BoolParameter.find(KEY)->second;    };


// The virtual getCauchyStress function that returns an error
int Material::getCauchyStress(  std::array<double, 6> &CauchyStress,
                                std::array<double, 9> &DeformationGradient,
                                std::map<std::string, double> &MaterialHistory,
                                std::map<std::string, int>    &IntegerMaterialIO,
                                std::map<std::string, double> &DoubleMaterialIO) {
                                std::cerr << "Error: Call of getCauchyStress() - function from Base Class ELSE::MPM::Material." << std::endl;
                                return -1;
                              };

}
}
