// 2D Material Point Method

#include <MPMProcess.hpp>
#include <MPMOutputVTK.hpp>
#include <MPM_Particle.hpp>
#include <MPM_GridNode.hpp>
//#include <MPM_OwnLib.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>

//declare function
void ReadParticleAreaData(double vec[]);
void ReadParticlePosition(double vec[][2]);
void TestVTUExport(double *xvecstart, int numx);


//------------------------------------------ MAIN ---------------------------------------------
int main()
{
    std::cout << "_______Welcome to MPM2D!______\n";
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now(); // Capture time at program start
    // First read in Particle Area data:
    // known particle len = 261
    //std::vector<double> Vp;
    double dens = 10;
    int i;
    // Read In The Problem:
    int nParticle = 216;          // Known number of particle
    //MPMParticle* GlobalParticleContainer = new MPMParticle[nParticle];
    std::vector<MPMParticle> GlobalParticleContainer;

    int NoGridNodes = 576;             // Number of GridNodes
    MPMGridNode* GlobalGridNodeContainer = new MPMGridNode[NoGridNodes];

    int nElements = 529;          // Number of GridElements
    double Vp[nParticle];         // Particle Volume Array
    double Xp[nParticle][2];      // Particle Position
    double XI[NoGridNodes][2];         // Particle Position
    //int nElements[nElements][4];  // Grid Element Connectivity

    ReadParticleAreaData(Vp);
    ReadParticlePosition(Xp);
    for(i=0;i<nParticle;i++){
      GlobalParticleContainer.push_back(MPMParticle(i+1,Xp[i][0],Xp[i][1],0.0,Vp[i],dens)); // put new born particle in the global container
    }
    // Hover over all particles in the global container and get a report
    std::vector<MPMParticle>::iterator v = GlobalParticleContainer.begin();
    while( v != GlobalParticleContainer.end()) {
    std::cout << "value of volume = " << (*v).Volume << std::endl;
    v++;
    }
    v = GlobalParticleContainer.begin();
    while( v != GlobalParticleContainer.end()) {
    (*v).Report();
    v++;
    }
    //for(i=0;i<nParticle;i++) std::cout << " Particle " << i+1 << " X=( " << Xp[i][0] << " , "<< Xp[i][1] << ")" << std::endl;

    // code for look at an array
    //for(i=0;i<261;i++) std::cout << i << "  " << Vp[i] << std::endl;

    MPMProcess MyProcess1; // Create a Process class (terminates automatically)
    MPMParticle MyFirstParticle(0.0,1.0,2.0);
    MPMOutputVTK MyOutput; // Create a MPMOutputVTK class (terminates automatically)
    //MyFirstParticle.Coordinate = [0.0, 0.0, 0.0];
    MyOutput.WriteVTK();
    MyFirstParticle.Mass = 10.0;
    //MyFirstParticle.Report();
    double x[4] = {0.0,1.0,1.0,0.0};
    double y[4] = {0.0,0.0,1.0,1.0};
    double z[4] = {0.0,0.0,0.0,0.0};
    TestVTUExport(x, 4);

    std::vector<double> xx;

    std::cout << "___________ The End __________\n";
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now(); // Capture time at program start
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout << "Runtime: " << duration*10e-6 << "s" << std::endl;
    //delete[] GlobalParticleContainer;
    return 0;
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

void ReadParticlePosition(double vec[][2]) {
  // create file object
  std::ifstream ip("/Users/sash/mpm_2d/data/two_discs_particle.txt"); // defiition of a ifstream object from the std space

  // check wether file could be loaded
  if(!ip.is_open()){
    std::cout << "ERROR: Cannot Open Particle Position File" << std::endl;
  } else {
    std::cout << "Start Reading Particle Position File" << std::endl;
  }

  //declare variables
  std::string ParticleX;
  std::string ParticleY;
  std::string::size_type nX;
  std::string::size_type nY;


  //loop over lines of file
  int i = 0;
  while(ip.good()){
    getline(ip, ParticleX, ',');
    getline(ip, ParticleY, '\n');
    try{
      vec[i][0] = std::stod(ParticleX, &nX);
      vec[i][1] = std::stod(ParticleY, &nY);
      i++;
    } catch (const std::exception& e) {
      std::cout << "End Reading Particle Position File" << std::endl;
    }
  }
}

void TestVTUExport(double *xvecstart, int numx){
  for(int m=0;m<numx;m++){
  std::cout << m+1 << "th X coordinate:" << *(xvecstart+m) << std::endl;
  }
}
