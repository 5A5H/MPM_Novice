#include <MPM_Particle.hpp>
#include <iostream>
#include <iomanip>

double MPMParticle::getMass(void){
  return Mass;
}

void MPMParticle::Report(void){
  std::cout << "------------------------------------\n";
  std::cout << "This is Particle :"  << std::setw(18) << ID << std::endl;
  std::cout << "Coordinate  X :"  << std::setw(21) << Coordinate[0] << std::endl;
  std::cout << "            Y :"  << std::setw(21) << Coordinate[1] << std::endl;
  std::cout << "            Z :"  << std::setw(21) << Coordinate[2] << std::endl;
  std::cout << "Mass          :"  << std::setw(21) << Mass << std::endl;
  std::cout << "Volume        :"  << std::setw(21) << Coordinate[0] << std::endl;
}
