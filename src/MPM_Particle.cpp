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
  for (int i=0;i<20;i++) h[i]=0.0;
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
  for (int i=0;i<20;i++) h[i]=0.0;
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
  for (int i=0;i<20;i++) h[i]=0.0;

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
  std::cout << "Sig11         :"  << std::setw(21) << Sig[0][0] << std::endl;
  std::cout << "Sig12         :"  << std::setw(21) << Sig[0][1] << std::endl;
  std::cout << "Sig13         :"  << std::setw(21) << Sig[0][2] << std::endl;
  std::cout << "Sig21         :"  << std::setw(21) << Sig[1][0] << std::endl;
  std::cout << "Sig22         :"  << std::setw(21) << Sig[1][1] << std::endl;
  std::cout << "Sig23         :"  << std::setw(21) << Sig[1][2] << std::endl;
  std::cout << "Sig31         :"  << std::setw(21) << Sig[2][0] << std::endl;
  std::cout << "Sig32         :"  << std::setw(21) << Sig[2][1] << std::endl;
  std::cout << "Sig33         :"  << std::setw(21) << Sig[2][2] << std::endl;
  std::cout << "F11         :"  << std::setw(21) << F[0][0] << std::endl;
  std::cout << "F12         :"  << std::setw(21) << F[0][1] << std::endl;
  std::cout << "F13         :"  << std::setw(21) << F[0][2] << std::endl;
  std::cout << "F21         :"  << std::setw(21) << F[1][0] << std::endl;
  std::cout << "F22         :"  << std::setw(21) << F[1][1] << std::endl;
  std::cout << "F23         :"  << std::setw(21) << F[1][2] << std::endl;
  std::cout << "F31         :"  << std::setw(21) << F[2][0] << std::endl;
  std::cout << "F32         :"  << std::setw(21) << F[2][1] << std::endl;
  std::cout << "F33         :"  << std::setw(21) << F[2][2] << std::endl;
  std::cout << "L11         :"  << std::setw(21) << L[0][0] << std::endl;
  std::cout << "L12         :"  << std::setw(21) << L[0][1] << std::endl;
  std::cout << "L13         :"  << std::setw(21) << L[0][2] << std::endl;
  std::cout << "L21         :"  << std::setw(21) << L[1][0] << std::endl;
  std::cout << "L22         :"  << std::setw(21) << L[1][1] << std::endl;
  std::cout << "L23         :"  << std::setw(21) << L[1][2] << std::endl;
  std::cout << "L31         :"  << std::setw(21) << L[2][0] << std::endl;
  std::cout << "L32         :"  << std::setw(21) << L[2][1] << std::endl;
  std::cout << "L33         :"  << std::setw(21) << L[2][2] << std::endl;
}

// Post-Processing Interface
// int: 0 -> key not defined !
// int: 1 -> scalar integer
// int: 2 -> scalar double
// int: 3 -> vector[3] double
// int: 4 -> tensor[9] double
int MPMParticle::GetPost(std::string KEY, std::array<double,9> &RealOut, std::array<int,9> &IntOut){
  if (KEY=="X"){RealOut[0]=X[0];RealOut[1]=X[1];RealOut[2]=X[2];return 3;}
  if (KEY=="V"){RealOut[0]=V[0];RealOut[1]=V[1];RealOut[2]=V[2];return 3;}
  if (KEY=="Mass"){RealOut[0]=Mass;return 2;}
  if (KEY=="J"){double detF;ELSE::Conti::Det(F,detF);RealOut[0]=detF;return 2;}
  if (KEY=="SigMises"){double sigmises;ELSE::Conti::VonMisesStress(Sig,sigmises);RealOut[0]=sigmises;return 2;}
  if (KEY=="TauMises"){RealOut[0]=MateData[3];return 2;}
  if (KEY=="MaterialState"){IntOut[0]= (MateData[0]<100)? 0:1;return 1;}
  if (KEY=="MaterialStatus"){IntOut[0]= (MateData[1]<100)? 0:1;return 1;}
  if (KEY=="MaterialIterations"){IntOut[0]= int( MateData[2] );return 1;}
  return 0;
}
