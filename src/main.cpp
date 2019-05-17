// 2D Material Point Method

#include <MPMProcess.hpp>
//#include <MPM_OwnLib.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

//declare function
void ReadParticleAreaData(double vec[]);
void PrintDoubleVector(std::vector<double> vec);


//------------------------------------------ MAIN ---------------------------------------------
int main()
{
    std::cout << "_______Welcome to MPM2D!______\n";
    // First read in Particle Area data:
    // known particle len = 261
    //std::vector<double> Vp;
    int i;
    double Vp[261];     // Particle Volume Array
    double Xp[261][2];  // Particle Position
    ReadParticleAreaData(Vp);
    for(i=0;i<261;i++) std::cout << Vp[i] << std::endl;

    MPMProcess MyProcess1; // Create a Process class (terminates automatically)

    return 0;
    std::cout << "_______ The End ______\n";
}

void ReadParticleAreaData(double vec[]) {
  // create file object
  std::ifstream ip("/Users/sash/mpm_2d/data/two_discs_particleArea.txt"); // defiition of a ifstream object from the std space

  // check wether file could be loaded
  if(!ip.is_open()){
    std::cout << "ERROR: Cannot Open File" << std::endl;
  } else {
    std::cout << "Start Reading File" << std::endl;
  }

  //declare variables
  std::string ParticleArea;
  std::string::size_type n;
  double WholeArea;


  //loop over lines of file
  WholeArea = 0.0;
  int i = 0;
  while(ip.good()){
    getline(ip, ParticleArea, '\n');
    //the read must be in a trail such that it fails when the expected pattern gets interrupted e.g. by a blank last line
    try{
      WholeArea += std::stod(ParticleArea, &n);
      vec[i] = std::stod(ParticleArea, &n);
      i++;
    } catch (const std::exception& e) {
      std::cout << "End Reading File" << std::endl;
    }
  }
  // Print for control the whole area
  std::cout << "WholeArea " << WholeArea << std::endl;
}

void PrintDoubleVector(std::vector<double> vec){
  std::vector<double>::iterator v = vec.begin();
  while( v != vec.end()) {
  std::cout << "value of v = " << *v << std::endl;
  v++;
  }
}
