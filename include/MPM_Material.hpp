#include <vector>
#include <iostream>
#include <AceMaterials.hpp>

class MPMMaterial {
  public:
      MPMMaterial(int id);
      ~MPMMaterial();
      int ID;                                                                   // Material ID (integer for now)
      // double X[3];                       // The GridNode Spatial Coordinate
      // double V[3];                       // The GridNode Velocity
      // double Density;                             // The GridNode Material Density
      // double Vol;                              // The GridNode Volume
      // double Mass;                                // The GridNode Mass
      // double Momentum[3];
      // double InternalForce[3];

      void SetMaterialParameter(double InputMaterialParameter);
      double *getStresses(double F[9]);                                       // for now only F but for future extended with optional args. e.g. time dependent
      void Report(void);                          // A Member Function to print out a report of this object

  private:
      std::vector<double> MaterialParameter;                                    // Contains arbitrary double material parameters which how to read is defined in the material subroutine of choice
};

MPMMaterial::~MPMMaterial(){}
MPMMaterial::MPMMaterial(int id){
  switch (id){ // ! muss wohl als loop umgestaltet werden und ein vairabler vector heallt dann die validen materialien
    case 1:
      std::cout << "Default Material Created (2d linear elastic) "  << std::endl;
      ID=id;
      break;

    case 2:
      std::cout << "AceGen material Created"  << std::endl;
      std::cout << " - linear elastic Hooke's law"  << std::endl;
      std::cout << " - plane strain assumption"  << std::endl;
      ID=id;
      break;

    case 3:
      std::cout << "AceGen material Created"  << std::endl;
      std::cout << " - linear elastic Hooke's law"  << std::endl;
      std::cout << " - plane stress assumption"  << std::endl;
      ID=id;
      break;

    default:
      std::cout << "Waring: Material not defined!"  << std::endl;
  }
}

void MPMMaterial::SetMaterialParameter(double InputMaterialParameter){
    MaterialParameter.push_back(InputMaterialParameter);
}

                                                                                //      F[0], F[1], F[2], F[3], F[4], F[5], F[6], F[7], F[8]
double *MPMMaterial::getStresses(double F[9]){                                  // F = [F11,  F12,  F13,  F21,  F22,  F23,  F31,  F32,  F33]
                                                                                //
  static double Sig[9] = {0,0,0,0,0,0,0,0,0};                                   // need to be defined as static to get out of scope of this function
  double Eps[6] = {0,0,0,0,0,0};
  double v[10];                                                                 // vector with auxillary variables (at least for default material)
  double vaux[100];                                                             // vector with auxillary variables (at least for acereturn)
  switch (ID){
    case 1:
          // Linear Elastic Default Material 2D
          if (MaterialParameter.size()!=2) {std::cout << "Warning: Material Parameter not set properly! "  << std::endl; break;}
          v[0] = MaterialParameter[0]; // Emod
          v[1] = MaterialParameter[1]; // nu
          v[2] = (v[0]*v[1]) / ( (1.0+v[1])*(1.0 - 2.0* v[1]) );  // lambda
          v[3] = v[0] / (2.0*(1.0 + v[1]));                       // mue
          Eps[0] = F[0]-1.0; // Eps11
          Eps[1] = F[4]-1.0; // Eps22
          Eps[2] = 0.5 * (F[1]+F[3]); // Eps12
          Sig[0] = (v[2]+2.0*v[3]) * Eps[0] + v[2] * Eps[1];  //Sig11
          Sig[4] = (v[2]+2.0*v[3]) * Eps[1] + v[2] * Eps[0];  //Sig22
          Sig[8] = v[2] * (Eps[0]+Eps[1]);                    //Sig33
          Sig[1] = 2.0*v[3] * Eps[2];                         //Sig12
          Sig[3] = Sig[1];                                    //Sig21

          break;

    case 2:
          // Call AceGen Material Subroutine: Linear Elastic Hookes Law with Plain Strain assumption
          if (MaterialParameter.size()!=2) {std::cout << "Warning: Material Parameter not set properly! "  << std::endl; break;}
          v[0] = MaterialParameter[0]; // Emod
          v[1] = MaterialParameter[1]; // nu
          Eps[0] = F[0]-1.0; // Eps11
          Eps[1] = F[4]-1.0; // Eps22
          Eps[2] = 0.5 * (F[1]+F[3]); // Eps12
          SmallStrainHookePlaneStrain2D(v, Eps, vaux);                          // call this way possible even if matrix dimensions dont agree, as eitherway only the pointer gets passed and the subroutine does not access beyond the dimensions of the arrays
          Sig[0] = vaux[0];
          Sig[4] = vaux[1];
          Sig[1] = vaux[3];
          Sig[3] = Sig[1];
          // note: this is a minmal version, a mature subroutine shall return the full stress tensor, but also gets the full deformation gradient !

          break;

    case 3:
          // Call AceGen Material Subroutine: Linear Elastic Hookes Law with Plain Strain assumption
          if (MaterialParameter.size()!=2) {std::cout << "Warning: Material Parameter not set properly! "  << std::endl; break;}
          v[0] = MaterialParameter[0]; // Emod
          v[1] = MaterialParameter[1]; // nu
          Eps[0] = F[0]-1.0; // Eps11
          Eps[1] = F[4]-1.0; // Eps22
          Eps[2] = 0.5 * (F[1]+F[3]); // Eps12
          SmallStrainHookePlaneStress2D(v, Eps, vaux);                          // call this way possible even if matrix dimensions dont agree, as eitherway only the pointer gets passed and the subroutine does not access beyond the dimensions of the arrays
          Sig[0] = vaux[0];
          Sig[4] = vaux[1];
          Sig[1] = vaux[3];
          Sig[3] = Sig[1];
          // note: this is a minmal version, a mature subroutine shall return the full stress tensor, but also gets the full deformation gradient !

          break;

    default:
    std::cout << "Waring: Material not defined!"  << std::endl;
  }                                                                             //        Sig[0], Sig[1], Sig[2], Sig[3], Sig[4], Sig[5], Sig[6], Sig[7], Sig[8]
  return Sig;                                                                   // Sig = [Sig11,  Sig12,  Sig13,  Sig21,  Sig22,  Sig23,  Sig31,  Sig32,  Sig33]
}

void MPMMaterial::Report(){
  std::cout << "---- Material Report ----" << std::endl;
  std::cout << "Material ID is: " << ID << std::endl;
  if (ID == 1 || ID == 2 || ID == 3){
    std::cout << "Young's modulus: " << MaterialParameter[0] << std::endl;
    std::cout << "Possion's ratio: " << MaterialParameter[1] << std::endl;
  }
  std::cout << "---- Material Report ----" << std::endl;
}
