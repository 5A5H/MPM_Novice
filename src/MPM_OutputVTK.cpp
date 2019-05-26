#include <MPM_OutputVTK.hpp>
#include <iostream>
#include <fstream>
#include <string>

// Constructor of the MPMProcess
MPMOutputVTK::MPMOutputVTK()
{
  std::cout << "VTK Started" << std::endl;
  FileName = "/Users/sash/mpm_2d/data/vpuout_1.vtp";  //For now hard wired
  vtpfile.open(FileName); //and associates the stream object to the file
}
// Destructor of the MPMProcess
MPMOutputVTK::~MPMOutputVTK()
{
  std::cout << "VTK Terminated" << std::endl;
}

int MPMOutputVTK::WriteVTK()
{
  vtpfile << "<VTKFile  type=\"PolyData\"  version=\"0.1\" >\n";
  vtpfile << "<PolyData>\n";
  vtpfile << "<Piece  NumberOfPoints=\"216\"  NumberOfVerts=\"0\"  NumberOfLines=\"0\"  NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
  vtpfile << "<Points>\n";
  WriteDataArray();
  vtpfile << "</Points>\n";
  vtpfile << "</Piece>\n";
  vtpfile << "</PolyData>\n";
  vtpfile << "</VTKFile>\n";
  vtpfile.close();
  return 0;
}

int MPMOutputVTK::WriteDataArray(){
  vtpfile << "<DataArray  type=\"Float64\"  NumberOfComponents=\"3\"  format=\"ascii\" >\n";
  // for(i=0;i<nParticle;i++){
  //   vtpfile << Xp[i][0] << "  " << Xp[i][1] << "  " << 0.0 << "\n";
  // }
  vtpfile << "</DataArray>\n";
  return 0;
}
