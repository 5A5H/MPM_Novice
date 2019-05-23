#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

void ReadParticle(std::string FileName, std::vector<MPMParticle> &ParticleContainer);
void ReadGridNodes(std::string FileName, std::vector<MPMGridNode> &GridNodeContainer);
void ReadGridElementsQ4(std::string FileName, std::vector<MPMGridElement> &GridElementContainer);

void ReadParticle(std::string FileName, std::vector<MPMParticle> &ParticleContainer) {
  std::ifstream INFile(FileName, std::ios::in);
  if (INFile.is_open()) {
    std::string line;
    int i = 0;
    double x,y,z,vol;
    while (std::getline(INFile, line)) {
      std::size_t pos = 0;
      // Find X-Coordinate
      pos = line.find(',');
      std::istringstream( line.substr(0,pos) ) >> x;
      //Find Y-Coordinae
      line.erase(0,pos+1);
      pos = line.find(',');
      std::istringstream( line.substr(0,pos) ) >> y;
      //Find Z-Coordinae
      line.erase(0,pos+1);
      pos = line.find(',');
      std::istringstream( line.substr(0,pos) ) >> z;
      //Find Volume
      line.erase(0,pos+1);
      pos = line.find(',');
      std::istringstream( line.substr(0,pos) ) >> vol;
      // Give birth to a particle
      ParticleContainer.push_back(   MPMParticle(i, x, y, z, vol)   );
      i++;
    }
    INFile.close();
  }
  else {
        std::cerr << "Unable to open file :"<< FileName << "\n";
  }
}

void ReadGridNodes(std::string FileName, std::vector<MPMGridNode> &GridNodeContainer) {
  std::ifstream INFile(FileName, std::ios::in);
  if (INFile.is_open()) {
    std::string line;
    int i = 0;
    double x,y;
    while (std::getline(INFile, line)) {
      std::size_t pos = 0;
      // Find X-Coordinate
      pos = line.find(',');
      std::istringstream( line.substr(0,pos) ) >> x;
      //Find Y-Coordinae
      line.erase(0,pos+1);
      pos = line.find(',');
      std::istringstream( line.substr(0,pos) ) >> y;
      // Give birth to a node
      GridNodeContainer.push_back(   MPMGridNode(i, x, y, 0.0)   );
      i++;
    }
    INFile.close();
  }
  else {
        std::cerr << "Unable to open file :"<< FileName << "\n";
  }
}

void ReadGridElementsQ4(std::string FileName, std::vector<MPMGridElement> &GridElementContainer) {
  std::ifstream INFile(FileName, std::ios::in);
  if (INFile.is_open()) {
    std::string line;
    int i = 0;
    int n1,n2,n3,n4;
    while (std::getline(INFile, line)) {
      std::size_t pos = 0;
      // Find X-Coordinate
      pos = line.find(',');
      std::istringstream( line.substr(0,pos) ) >> n1;
      //Find Y-Coordinae
      line.erase(0,pos+1);
      pos = line.find(',');
      std::istringstream( line.substr(0,pos) ) >> n2;
      //Find Z-Coordinae
      line.erase(0,pos+1);
      pos = line.find(',');
      std::istringstream( line.substr(0,pos) ) >> n3;
      //Find Z-Coordinae
      line.erase(0,pos+1);
      pos = line.find(',');
      std::istringstream( line.substr(0,pos) ) >> n4;
      // Give birth to a node
      GridElementContainer.push_back(   MPMGridElement(i, n1, n2, n3, n4)   );
      i++;
    }
    INFile.close();
  }
  else {
    std::cerr << "Unable to open file :"<< FileName << "\n";
  }
}
