#include <MPM_Particle.hpp>
#include <iostream>
#include <iomanip>

// Constructor of the Particle class
MPMParticle::MPMParticle(double x, double y, double z){
  ID = 0;
  Coordinate[0] = x;
  Coordinate[1] = y;
  Coordinate[2] = z;
  Density = 1;
  Volume = 1;
  Mass = Density*Volume;
}
MPMParticle::MPMParticle(){
  ID = 0;
  Density = 1;
  Volume = 1;
  Mass = Density*Volume;
}
MPMParticle::MPMParticle(int id, double x, double y, double z, double vol, double dens){
  ID = id;
  Coordinate[0] = x;
  Coordinate[1] = y;
  Coordinate[2] = z;
  Density = dens;
  Volume = vol;
  Mass = Density*Volume;
}

// Destructor of the Particle class
MPMParticle::~MPMParticle(){
}

// A member function to return the mass (all dummy)
double MPMParticle::getMass(void){
  return Mass;
}

// particle data report to prompt
void MPMParticle::Report(void){
  std::cout << "------------------------------------\n";
  std::cout << "This is Particle :"  << std::setw(18) << ID << std::endl;
  std::cout << "Coordinate  X :"  << std::setw(21) << Coordinate[0] << std::endl;
  std::cout << "            Y :"  << std::setw(21) << Coordinate[1] << std::endl;
  std::cout << "            Z :"  << std::setw(21) << Coordinate[2] << std::endl;
  std::cout << "Mass          :"  << std::setw(21) << Mass << std::endl;
  std::cout << "Volume        :"  << std::setw(21) << Volume << std::endl;
}
