#include <MPM_GridNode.hpp>
#include <iostream>
#include <iomanip>

// Constructor of the Particle class
MPMGridNode::MPMGridNode(int id, double x, double y, double z){
  ID = id;
  X[0] = x;
  X[1] = y;
  X[2] = z;
  Density = 0;
  Vol = 0;
  Mass = 0;
  Momentum[0] = 0.0;
  Momentum[1] = 0.0;
  Momentum[2] = 0.0;
  InternalForce[0] = 0;
  InternalForce[1] = 0;
  InternalForce[2] = 0;
  Force[0] = 0;
  Force[1] = 0;
  Force[2] = 0;
  V[0] = 0;
  V[1] = 0;
  V[2] = 0;
}
MPMGridNode::MPMGridNode(){
  ID = -1;
  X[0] = 0.0;
  X[1] = 0.0;
  X[2] = 0.0;
  Density = 0;
  Vol = 0;
  Mass = 0;
  Momentum[0] = 0.0;
  Momentum[1] = 0.0;
  Momentum[2] = 0.0;
  InternalForce[0] = 0;
  InternalForce[1] = 0;
  InternalForce[2] = 0;
  V[0] = 0;
  V[1] = 0;
  V[2] = 0;
  Force[0] = 0;
  Force[1] = 0;
  Force[2] = 0;
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
  // Reset Grid Forces
  Force[0] = 0;
  Force[1] = 0;
  Force[2] = 0;
}
