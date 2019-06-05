#ifndef _MPM_GRID_NODE_
#define _MPM_GRID_NODE_

class MPMGridNode {
  public:
      MPMGridNode();                              // The GridNode Class Constructor
      MPMGridNode(int id, double x, double y, double z);  // The GridNode Class Constructor requires coordinate
      ~MPMGridNode();                             // The GridNode Class Destructor
      int ID;                                     // The GridNode Id
      double X[3];                       // The GridNode Spatial Coordinate
      double V[3];                       // The GridNode Velocity
      double Density;                             // The GridNode Material Density
      double Vol;                              // The GridNode Volume
      double Mass;                                // The GridNode Mass
      double Momentum[3];
      double InternalForce[3];
      double ExternalForce[3];
      double Force[3];

      void Reset();
      void Report(void);                          // A Member Function to print out a report of this object

  private:
};

#endif
