// 2D Material Point Method

#include <MPMProcess.hpp>
#include <MPMOutputVTK.hpp>
#include <MPM_Particle.hpp>
//#include <MPM_OwnLib.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

//declare function
void ReadParticleAreaData(double vec[]);
void ReadParticlePosition(double vec[][2]);


//------------------------------------------ MAIN ---------------------------------------------
int main()
{
    std::cout << "_______Welcome to MPM2D!______\n";
    // First read in Particle Area data:
    // known particle len = 261
    //std::vector<double> Vp;
    int i;
    int nParticle = 216;        //Known number of particle
    double Vp[nParticle];       // Particle Volume Array
    double Xp[nParticle][2];    // Particle Position
    ReadParticleAreaData(Vp);
    ReadParticlePosition(Xp);
    //for(i=0;i<nParticle;i++) std::cout << " Particle " << i+1 << " X=( " << Xp[i][0] << " , "<< Xp[i][1] << ")" << std::endl;

    // code for look at an array
    //for(i=0;i<261;i++) std::cout << i << "  " << Vp[i] << std::endl;

    // write to vtu file
    std::ofstream vtpfile;
    vtpfile.open("/Users/sash/mpm_2d/data/vpuout_1.vtp");
    vtpfile << "<VTKFile  type=\"PolyData\"  version=\"0.1\" >\n";
    vtpfile << "<PolyData>\n";
    vtpfile << "<Piece  NumberOfPoints=\"216\"  NumberOfVerts=\"0\"  NumberOfLines=\"0\"  NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
    vtpfile << "<Points>\n";
    vtpfile << "<DataArray  type=\"Float64\"  NumberOfComponents=\"3\"  format=\"ascii\" >\n";
    for(i=0;i<nParticle;i++){
      vtpfile << Xp[i][0] << "  " << Xp[i][1] << "  " << 0.0 << "\n";
    }
    vtpfile << "</DataArray>\n";
    vtpfile << "</Points>\n";
    vtpfile << "</Piece>\n";
    vtpfile << "</PolyData>\n";
    vtpfile << "</VTKFile>\n";
    vtpfile.close();


    MPMProcess MyProcess1; // Create a Process class (terminates automatically)
    MPMOutputVTK MyOutput; // Create a MPMOutputVTK class (terminates automatically)
    MPMParticle MyFirstParticle;
    //MyFirstParticle.Coordinate = [0.0, 0.0, 0.0];
    MyFirstParticle.Mass = 10.0;
    std::cout << "MyFirstParticle has Mass: " << MyFirstParticle.getMass() << std::endl;
    MyFirstParticle.Report();

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
