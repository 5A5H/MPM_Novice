#include <MPM_Particle.hpp>
#include <iostream>
#include <iomanip>

// Constructor of the Particle class
MPMParticle::MPMParticle(double x, double y, double z){
  ID = 0;
  X[0] = x;
  X[1] = y;
  X[2] = z;
  Vol = 1;
  Mass = 0;
  V[0] = 0.0;
  V[1] = 0.0;
  V[2] = 0.0;
  Stress[0] = 0.0;
  Stress[1] = 0.0;
  Stress[2] = 0.0;
  Deformation[0] = 1.0;
  Deformation[1] = 0.0;
  Deformation[2] = 0.0;
  Deformation[3] = 1.0;
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      Sig[i][j] = 0e0;
      F[i][j] = 0e0;
      L[i][j] = 0e0;
    }
    F[i][i] = 1e0;
    b[i] = 0e0;
  }
}
MPMParticle::MPMParticle(){
  ID = 0;
  Vol = 1;
  Mass = 0;
  V[0] = 0.0;
  V[1] = 0.0;
  V[2] = 0.0;
  Stress[0] = 0.0;
  Stress[1] = 0.0;
  Stress[2] = 0.0;
  Deformation[0] = 1.0;
  Deformation[1] = 0.0;
  Deformation[2] = 0.0;
  Deformation[3] = 1.0;
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      Sig[i][j] = 0e0;
      F[i][j] = 0e0;
      L[i][j] = 0e0;
    }
    F[i][i] = 1e0;
    b[i] = 0e0;
  }
}
MPMParticle::MPMParticle(int id, double x, double y, double z, double vol){
  ID = id;
  X[0] = x;
  X[1] = y;
  X[2] = z;
  Vol = vol;
  Mass = 0.0;
  V[0] = 0.0;
  V[1] = 0.0;
  V[2] = 0.0;
  Stress[0] = 0.0;
  Stress[1] = 0.0;
  Stress[2] = 0.0;
  Deformation[0] = 1.0;
  Deformation[1] = 0.0;
  Deformation[2] = 0.0;
  Deformation[3] = 1.0;

  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      Sig[i][j] = 0e0;
      F[i][j] = 0e0;
      L[i][j] = 0e0;
    }
    F[i][i] = 1e0;
    b[i] = 0e0;
  }

}

// Destructor of the Particle class
MPMParticle::~MPMParticle(){
}

// A member function to return the mass (all dummy)
double MPMParticle::getMass(void){
  return Mass;
}

// Check for check
bool MPMParticle::checkNAN(){
  if (X[0]!=X[0] || X[1]!=X[1] || X[2]!=X[2]){
    std::cout << std::endl;
    std::cerr << "Unfortunately NAN detected in Particle :" << ID << "\n";
    return true;
  } else {
    return false;
 }
}

// particle data report to prompt
void MPMParticle::Report(void){
  std::cout << "------------------------------------\n";
  std::cout << "This is Particle :"  << std::setw(18) << ID << std::endl;
  std::cout << "Coordinate  X :"  << std::setw(21) << X[0] << std::endl;
  std::cout << "            Y :"  << std::setw(21) << X[1] << std::endl;
  std::cout << "            Z :"  << std::setw(21) << X[2] << std::endl;
  std::cout << "Mass          :"  << std::setw(21) << Mass << std::endl;
  std::cout << "Volume        :"  << std::setw(21) << Vol << std::endl;
}
