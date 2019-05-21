class MPMParticle {
  public:
      MPMParticle();                              // The Particle Class Constructor
      MPMParticle(double x, double y, double z);  // The Particle Class Constructor requires coordinate
      MPMParticle(int id, double x, double y, double z, double vol, double dens);
      ~MPMParticle();                             // The Particle Class Destructor
      int ID;                                     // The Particles Id
      double X[3];                       // The Particles Spatial Coordinate
      double V[3];                         // The Particles Velocity
      double Density;                             // The Particles Material Density
      double Vol;                              // The Particles Volume
      double Mass;                                // The Particles Mass
      double Stress;                                // The Particles Stress {Sig11,Sig22,Sig12}
      double Deformation;                                // The Particles Deformation {F11,F12,F21,F22}

      double getMass(void);                       // A Member Function that returns the Current Mass of the Particle  MemberFunction: hass access to all data of the object
      void Report(void);                          // A Member Function to print out a report of this object

  private:
};
