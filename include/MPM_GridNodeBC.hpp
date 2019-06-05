#include <iostream>
#include <vector>

class MPMGridNodeBC {
  public:
      MPMGridNodeBC();                              // The GridNode Class Constructor
      ~MPMGridNodeBC();                             // The GridNode Class Destructor

      // int ID;                                     // The GridNode Id
      // double X[3];                       // The GridNode Spatial Coordinate
      // double V[3];                       // The GridNode Velocity
      // double Density;                             // The GridNode Material Density
      // double Vol;                              // The GridNode Volume
      // double Mass;                                // The GridNode Mass
      // double Momentum[3];
      // double InternalForce[3];
      //
      // void Reset();
      // void Report(void);                          // A Member Function to print out a report of this object
      void setBC(std::string type, std::string affectedProperty, double prescribedValue);
      void addGridNode(int id);
      void applyBC(std::vector<MPMGridNode> &AllNodeContainer);

  private:
    std::vector<int> AffectedNodes;                                             // vector of node ids where the bc has effect
};

MPMGridNodeBC::MPMGridNodeBC(){}
MPMGridNodeBC::~MPMGridNodeBC(){}

void MPMGridNodeBC::setBC(std::string type, std::string affectedProperty, double prescribedValue){
  bool DetailedOutput = true;
  if(type=="EssentialBC"){
    if (DetailedOutput) std::cout << "Essential BC :  " << affectedProperty << " -> " << prescribedValue << "  created." << std::endl ;
  }
}

void MPMGridNodeBC::addGridNode(int id){
  AffectedNodes.push_back(id);
}

void MPMGridNodeBC::applyBC(std::vector<MPMGridNode> &AllNodeContainer){
  //std::cout << AffectedNodes.size() << std::endl;
  int j;
  for (j=0;j<AffectedNodes.size();j++){
     //AllNodeContainer[AffectedNodes[j]].V[0] = 0.0;
     AllNodeContainer[AffectedNodes[j]].V[1] = 0.0;
     //AllNodeContainer[AffectedNodes[j]].V[2] = 0.0;
     //AllNodeContainer[AffectedNodes[j]].Momentum[0] = 0.0;
     AllNodeContainer[AffectedNodes[j]].Momentum[1] = 0.0;
     //AllNodeContainer[AffectedNodes[j]].Momentum[2] = 0.0;
     //AllNodeContainer[AffectedNodes[j]].InternalForce[0] = 0.0;
     AllNodeContainer[AffectedNodes[j]].InternalForce[1] = 0.0;
     //AllNodeContainer[AffectedNodes[j]].InternalForce[2] = 0.0;
     //std::cout << "set bc on node: " << AffectedNodes[j] << "   X: " << AllNodeContainer[AffectedNodes[j]].X[0] << ", " << AllNodeContainer[AffectedNodes[j]].X[1] << ", " << AllNodeContainer[AffectedNodes[j]].X[2] << std::endl;
   }
}
