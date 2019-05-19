#include <MPM_GridNode.hpp>
#include <iostream>
#include <iomanip>

// Constructor of the Particle class
MPMGridNode::MPMGridNode(double x, double y, double z){
  ID = 0;
  Coordinate[0] = x;
  Coordinate[1] = y;
  Coordinate[2] = z;
  Density = 1;
  Volume = 1;
  Mass = Density*Volume;
}
MPMGridNode::MPMGridNode(){
  ID = 0;
  Density = 1;
  Volume = 1;
  Mass = Density*Volume;
}

// Destructor of the Particle class
MPMGridNode::~MPMGridNode(){
}


// particle data report to prompt
void MPMGridNode::Report(void){
  std::cout << "------------------------------------\n";
  std::cout << "This is GridNode :"  << std::setw(18) << ID << std::endl;
  std::cout << "Coordinate  X :"  << std::setw(21) << Coordinate[0] << std::endl;
  std::cout << "            Y :"  << std::setw(21) << Coordinate[1] << std::endl;
  std::cout << "            Z :"  << std::setw(21) << Coordinate[2] << std::endl;
  std::cout << "Mass          :"  << std::setw(21) << Mass << std::endl;
  std::cout << "Volume        :"  << std::setw(21) << Volume << std::endl;
}
