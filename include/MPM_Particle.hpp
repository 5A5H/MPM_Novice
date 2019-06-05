#ifndef _MPM_PARTICLE_
#define _MPM_PARTICLE_


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

      bool checkNAN();
      double getMass(void);                       // A Member Function that returns the Current Mass of the Particle  MemberFunction: hass access to all data of the object
      void Report(void);                          // A Member Function to print out a report of this object

  private:
};

#endif
