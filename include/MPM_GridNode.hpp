class MPMGridNode {
  public:
      MPMGridNode();                              // The Particle Class Constructor
      MPMGridNode(double x, double y, double z);  // The Particle Class Constructor requires coordinate
      ~MPMGridNode();                             // The Particle Class Destructor
      int ID;                                     // The Particles Id
      double Coordinate[3];                       // The Particles Spatial Coordinate
      double Density;                             // The Particles Material Density
      double Volume;                              // The Particles Volume
      double Mass;                                // The Particles Mass

      void Report(void);                          // A Member Function to print out a report of this object

  private:
};
