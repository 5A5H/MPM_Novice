// 2D Material Point Method

#include <MPMProcess.hpp>
#include <MPMOutputVTK.hpp>
#include <MPM_Particle.hpp>
#include <MPM_GridNode.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>

//declare function
void ReadParticleAreaData(double vec[]);
void ReadParticlePosition(double vec[][2]);
void TestVTUExport(std::vector<MPMParticle> &OutParticleContainer);


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
    //std::cout << "value of volume = " << (*v).Volume << std::endl;
    v++;
    }
    v = GlobalParticleContainer.begin();
    while( v != GlobalParticleContainer.end()) {
    //(*v).Report();
    v++;
    }

    double sumvol = 0;
    v = GlobalParticleContainer.begin();
    while( v != GlobalParticleContainer.end()) {
    sumvol += (*v).Volume;
    v++;
    }
    std::cout << "Complete Particle Volume :" << sumvol << std::endl;

    // code for look at an array
    //for(i=0;i<261;i++) std::cout << i << "  " << Vp[i] << std::endl;

    MPMProcess MyProcess1; // Create a Process class (terminates automatically)
    MPMOutputVTK MyOutput; // Create a MPMOutputVTK class (terminates automatically)
    MyOutput.WriteVTK();
    TestVTUExport(GlobalParticleContainer);

    std::vector<double> xx;

    std::cout << "___________ The End __________\n";
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now(); // Capture time at program start
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout << "Runtime: " << duration*10e-6 << "s" << std::endl;
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

void TestVTUExport(std::vector<MPMParticle> &OutParticleContainer){                     //GlobalParticleContainerLength->GPCL #, int GPCL
  std::cout << "container len:" << OutParticleContainer.size() << std::endl;

  //collecting information
  //piece information
  int NumberOfPoints  =   OutParticleContainer.size();
  int NumberOfCells   =   0;


  //open a file stream for output
  std::ofstream OutputFile;
  OutputFile.open("/Users/sash/mpm_2d/data/vpuout_2.vtu", std::ios::out);
  // Write headder
  OutputFile << "<?xml version=\"1.0\" ?>" << std::endl;
  OutputFile << "<VTKFile byte_order=\"LittleEndian\" type=\"UnstructuredGrid\" version=\"0.1\">" << std::endl;
  OutputFile << "<UnstructuredGrid>" << std::endl;
  // Write Piece Headder
  OutputFile << "<Piece NumberOfCells=\"" << NumberOfCells << "\" NumberOfPoints=\"" << NumberOfPoints << "\">" << std::endl;
  // Write Points
  OutputFile << "<Points>" << std::endl;

  OutputFile << "<DataArray NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">";
  for(int i = 0; i < NumberOfPoints; i++){
    OutputFile << "  " << OutParticleContainer[i].Coordinate[0] ;
    OutputFile << "  " << OutParticleContainer[i].Coordinate[1] ;
    OutputFile << "  " << OutParticleContainer[i].Coordinate[2] ;
  }
  OutputFile << "</DataArray>" << std::endl;

  OutputFile << "</Points>" << std::endl;
  // Write Cells
  OutputFile << "<Cells>" << std::endl;

  OutputFile << "<DataArray Name=\"connectivity\" format=\"ascii\" type=\"Int32\">0</DataArray>" << std::endl;
  OutputFile << "<DataArray Name=\"offsets\" format=\"ascii\" type=\"Int32\">0</DataArray>" << std::endl;
  OutputFile << "<DataArray Name=\"types\" format=\"ascii\" type=\"UInt8\">1</DataArray>" << std::endl;

  OutputFile << "</Cells>" << std::endl;
  // Write Point Data
  OutputFile << "<PointData/>" << std::endl;
  // Write Cell Data
  OutputFile << "<CellData/>" << std::endl;
  // Write Piece foot
  OutputFile << "</Piece>" << std::endl;
  // Write foot
  OutputFile << "</UnstructuredGrid>" << std::endl;
  OutputFile << "</VTKFile>" << std::endl;
  //close the file stram
  OutputFile.close();
}


// void TestVTUExport(double *xvecstart, int numx){
//   for(int m=0;m<numx;m++){
//   std::cout << m+1 << "th X coordinate:" << *(xvecstart+m) << std::endl;
//   }
// }
