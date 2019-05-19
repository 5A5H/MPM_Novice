class MPMParticle {
  public:
      MPMParticle();                              // The Particle Class Constructor
      MPMParticle(double x, double y, double z);  // The Particle Class Constructor requires coordinate
      MPMParticle(int id, double x, double y, double z, double vol, double dens);
      ~MPMParticle();                             // The Particle Class Destructor
      int ID;                                     // The Particles Id
      double Coordinate[3];                       // The Particles Spatial Coordinate
      double Density;                             // The Particles Material Density
      double Volume;                              // The Particles Volume
      double Mass;                                // The Particles Mass

      double getMass(void);                       // A Member Function that returns the Current Mass of the Particle  MemberFunction: hass access to all data of the object
      void Report(void);                          // A Member Function to print out a report of this object

  private:
};
