#include <MPM_GridNode.hpp>
#include <iostream>
#include <iomanip>

// Constructor of the Particle class
MPMGridNode::MPMGridNode(double x, double y, double z){
  ID = 0;
  X[0] = x;
  X[1] = y;
  X[2] = z;
  Density = 1;
  Vol = 1;
  Mass = 0;
}
MPMGridNode::MPMGridNode(){
  ID = 0;
  Density = 1;
  Vol= 1;
  Mass = 0;
}

// Destructor of the Particle class
MPMGridNode::~MPMGridNode(){
}


// particle data report to prompt
void MPMGridNode::Report(void){
  std::cout << "------------------------------------\n";
  std::cout << "This is GridNode :"  << std::setw(18) << ID << std::endl;
  std::cout << "Coordinate  X :"  << std::setw(21) << X[0] << std::endl;
  std::cout << "            Y :"  << std::setw(21) << X[1] << std::endl;
  std::cout << "            Z :"  << std::setw(21) << X[2] << std::endl;
  std::cout << "Mass          :"  << std::setw(21) << Mass << std::endl;
  std::cout << "Volume        :"  << std::setw(21) << Vol << std::endl;
}

void MPMGridNode::Reset(){
  // Reset Nodal Mass
  Mass = 0;
  // Reset Nodal Momentum
  Momentum[0] = 0;
  Momentum[1] = 0;
  Momentum[2] = 0;
  // Reset Nodal Internal Forces
  InternalForce[0] = 0;
  InternalForce[1] = 0;
  InternalForce[2] = 0;
  // Reset Nodal Velocity
  V[0] = 0;
  V[1] = 0;
  V[2] = 0;
}
