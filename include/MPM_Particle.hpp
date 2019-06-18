#ifndef _MPM_PARTICLE_
#define _MPM_PARTICLE_

#include <math.h>
#include <string>
#include <vector>
#include <array>
#include <ELSE_ContiMechLibrary.hpp>

class MPMParticle {
  public:
      MPMParticle();                              // The Particle Class Constructor
      MPMParticle(double x, double y, double z);  // The Particle Class Constructor requires coordinate
      MPMParticle(int id, double x, double y, double z, double vol);
      ~MPMParticle();                             // The Particle Class Destructor




      int ID;                                     // The Particles Id
      double X[3];                       // The Particles Spatial Coordinate
      double V[3];                         // The Particles Velocity
      double Vol;                              // The Particles Volume
      double Mass;                                // The Particles Mass
      double Density;                              // The Particles Density
      double Stress[3];                                // The Particles Stress {Sig11,Sig22,Sig12}
      double Deformation[4];                                // The Particles Deformation {F11,F12,F21,F22}
      double Sig[3][3]; // cauchy stresses
      double F[3][3]; //deformation gradient
      double L[3][3]; //velocity gradient
      double b[3]; //body force
      int Elmt; //current element of the particle
      double h[20]; // material point history
      double MateData[20];//additional data form material


      bool checkNAN();
      double getMass(void);                       // A Member Function that returns the Current Mass of the Particle  MemberFunction: hass access to all data of the object
      void Report(void);                          // A Member Function to print out a report of this object
      double SigMises(){return sqrt(0.15e1*(2e0*Sig[0][1]*Sig[1][0]+2e0*Sig[0][2]*Sig[2][0]+2e0*Sig[1][2]*Sig[2][1]+pow(Sig[0][0]+((-Sig[0][0]-Sig[1][1]-Sig[2][2])/3e0),2)+pow(Sig[1][1]+((-Sig[0][0]-Sig[1][1]-Sig[2][2])/3e0),2)+pow(Sig[2][2]+((-Sig[0][0]-Sig[1][1]-Sig[2][2])/3e0),2)));}

      // Implementation of Post-Processing
      // Function that just return smth.
      // Input is a key
      // each function returns an integer for its information and modifies up to 9->trensor double values
      int GetPost(std::string KEY, std::array<double,9> &RealOut, std::array<int,9> &IntOut);
      // int: 0 -> key not defined !
      // int: 1 -> scalar integer
      // int: 2 -> scalar double
      // int: 3 -> vector[3] double
      // int: 4 -> tensor[9] double
  private:
};

//MateData
// 1->MaterialState
// 2->material convergence
// 3->no material iterations
// 4-> final taumises


#endif
