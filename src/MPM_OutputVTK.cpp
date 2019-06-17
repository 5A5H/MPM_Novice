#include <MPM_OutputVTK.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

MPMOutputVTK::MPMOutputVTK(){}

MPMOutputVTK::~MPMOutputVTK(){
  WritePVD();
}

void MPMOutputVTK::TestVTUParticleExport(std::string FileName, std::vector<MPMParticle> &OutParticleContainer){


  //collecting information
  //piece information
  int NumberOfParticles  = OutParticleContainer.size();


  //open a file stream for output
  std::ofstream OutputFile;
  OutputFile.open(FileName, std::ios::out);
  // Write headder
  OutputFile << "<?xml version=\"1.0\" ?>" << std::endl;
  OutputFile << "<VTKFile byte_order=\"LittleEndian\" type=\"UnstructuredGrid\" version=\"0.1\">" << std::endl;
  OutputFile << "<UnstructuredGrid>" << std::endl;

  // Piece 1 -> Particle
  // Write Piece Headder
  OutputFile << "<Piece NumberOfCells=\"" << NumberOfParticles << "\" NumberOfPoints=\"" << NumberOfParticles << "\">" << std::endl;
  // Write Points
  OutputFile << "<Points>" << std::endl;

  OutputFile << "<DataArray NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">";
  for(int i = 0; i < NumberOfParticles; i++){
    OutputFile << "  " << CuttOff(OutParticleContainer[i].X[0]) ;
    OutputFile << "  " << CuttOff(OutParticleContainer[i].X[1]) ;
    OutputFile << "  " << CuttOff(OutParticleContainer[i].X[2]) ;
  }
  OutputFile << "</DataArray>" << std::endl;

  OutputFile << "</Points>" << std::endl;
  // Write Cells
  OutputFile << "<Cells>" << std::endl;

  OutputFile << "<DataArray Name=\"connectivity\" format=\"ascii\" type=\"Int32\">";
  for(int i = 0; i < NumberOfParticles; i++){
    OutputFile << "  " << i ;
  }
  OutputFile << "</DataArray>" << std::endl;

  OutputFile << "<DataArray Name=\"offsets\" format=\"ascii\" type=\"Int32\">";
  for(int i = 0; i < NumberOfParticles; i++) OutputFile << "  " << (i+1)*1;
  OutputFile << "</DataArray>" << std::endl;


  OutputFile << "<DataArray Name=\"types\" format=\"ascii\" type=\"UInt8\">";
  for(int i = 0; i < NumberOfParticles; i++) OutputFile << "  " << 1;
  OutputFile << "</DataArray>" << std::endl;

  //OutputFile << "<DataArray Name=\"connectivity\" format=\"ascii\" type=\"Int32\">0</DataArray>" << std::endl;
  //OutputFile << "<DataArray Name=\"offsets\" format=\"ascii\" type=\"Int32\">0</DataArray>" << std::endl;
  //OutputFile << "<DataArray Name=\"types\" format=\"ascii\" type=\"UInt8\">1</DataArray>" << std::endl;

  OutputFile << "</Cells>" << std::endl;
  // Write Point Data
  OutputFile << "<PointData>" << std::endl;
    // Write Particle Volume
    OutputFile << "<DataArray Name=\"Volume\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
    for(int i = 0; i < NumberOfParticles; i++){
      OutputFile << "  " << CuttOff(OutParticleContainer[i].Vol) ;
    }
    OutputFile << "</DataArray>" << std::endl;
    // Write Particle Mass
    OutputFile << "<DataArray Name=\"Mass\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
    for (auto &Particle : OutParticleContainer) {
      OutputFile << "  " << CuttOff(Particle.Mass);
    }
    OutputFile << "</DataArray>" << std::endl;
/*
    // Write Particle Stresses
    OutputFile << "<DataArray Name=\"Sig11\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
    for (auto &Particle : OutParticleContainer) {
      OutputFile << "  " << CuttOff(Particle.Sig[0][0]);
    }
    OutputFile << "</DataArray>" << std::endl;
    OutputFile << "<DataArray Name=\"Sig12\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
    for (auto &Particle : OutParticleContainer) {
      OutputFile << "  " << CuttOff(Particle.Sig[0][1]);
    }
    OutputFile << "</DataArray>" << std::endl;
    OutputFile << "<DataArray Name=\"Sig13\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
    for (auto &Particle : OutParticleContainer) {
      OutputFile << "  " << CuttOff(Particle.Sig[0][2]);
    }
    OutputFile << "</DataArray>" << std::endl;
    OutputFile << "<DataArray Name=\"Sig21\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
    for (auto &Particle : OutParticleContainer) {
      OutputFile << "  " << CuttOff(Particle.Sig[1][0]);
    }
    OutputFile << "</DataArray>" << std::endl;
    OutputFile << "<DataArray Name=\"Sig22\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
    for (auto &Particle : OutParticleContainer) {
      OutputFile << "  " << CuttOff(Particle.Sig[1][1]);
    }
    OutputFile << "</DataArray>" << std::endl;
    OutputFile << "<DataArray Name=\"Sig23\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
    for (auto &Particle : OutParticleContainer) {
      OutputFile << "  " << CuttOff(Particle.Sig[1][2]);
    }
    OutputFile << "</DataArray>" << std::endl;
    OutputFile << "<DataArray Name=\"Sig31\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
    for (auto &Particle : OutParticleContainer) {
      OutputFile << "  " << CuttOff(Particle.Sig[2][0]);
    }
    OutputFile << "</DataArray>" << std::endl;
    OutputFile << "<DataArray Name=\"Sig32\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
    for (auto &Particle : OutParticleContainer) {
      OutputFile << "  " << CuttOff(Particle.Sig[2][1]);
    }
    OutputFile << "</DataArray>" << std::endl;
    OutputFile << "<DataArray Name=\"Sig33\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
    for (auto &Particle : OutParticleContainer) {
      OutputFile << "  " << CuttOff(Particle.Sig[2][2]);
    }
    OutputFile << "</DataArray>" << std::endl;
*/
    // Write VonMises Stresses
    OutputFile << "<DataArray Name=\"SigMises\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
    for (auto &Particle : OutParticleContainer) {
      double sigmises = Particle.SigMises();
      OutputFile << "  " << CuttOff(sigmises,1e10,-1);
    }
    OutputFile << "</DataArray>" << std::endl;

    // Write plastic indicator
    OutputFile << "<DataArray Name=\"alpha\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
    for (auto &Particle : OutParticleContainer) {
      OutputFile << "  " << CuttOff(Particle.h[9]);
    }
    OutputFile << "</DataArray>" << std::endl;



    // Write Particle Velocity
    OutputFile << "<DataArray Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">";
    for(int i = 0; i < NumberOfParticles; i++){
      OutputFile << "  " << CuttOff(OutParticleContainer[i].V[0]);
      OutputFile << "  " << CuttOff(OutParticleContainer[i].V[1]);
      OutputFile << "  " << CuttOff(OutParticleContainer[i].V[2]);
    }
    OutputFile << "</DataArray>" << std::endl;
  OutputFile << "</PointData>" << std::endl;
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

void MPMOutputVTK::TestVTUGridExport(
  std::string FileName,
  std::vector<MPMGridNode> &OutNodeContainer,
  std::vector<MPMGridElement> &OutElementContainer){
    {


      //collecting information
      //piece information
      int NumberOfNodes   =   OutNodeContainer.size();
      int NumberOfCells   =   OutElementContainer.size();


      //open a file stream for output
      std::ofstream OutputFile;
      OutputFile.open(FileName, std::ios::out);
      // Write headder
      OutputFile << "<?xml version=\"1.0\" ?>" << std::endl;
      OutputFile << "<VTKFile byte_order=\"LittleEndian\" type=\"UnstructuredGrid\" version=\"0.1\">" << std::endl;
      OutputFile << "<UnstructuredGrid>" << std::endl;

      // Piece 1 -> Particle
      // Write Piece Headder
      OutputFile << "<Piece NumberOfCells=\"" << NumberOfCells << "\" NumberOfPoints=\"" << NumberOfNodes << "\">" << std::endl;
      // Write Points
      OutputFile << "<Points>" << std::endl;

      OutputFile << "<DataArray NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">";
      for(int i = 0; i < NumberOfNodes; i++){
        OutputFile << "  " << OutNodeContainer[i].X[0] ;
        OutputFile << "  " << OutNodeContainer[i].X[1] ;
        OutputFile << "  " << OutNodeContainer[i].X[2] ;
      }
      OutputFile << "</DataArray>" << std::endl;

      OutputFile << "</Points>" << std::endl;
      // Write Cells
      OutputFile << "<Cells>" << std::endl;

          OutputFile << "<DataArray Name=\"connectivity\" format=\"ascii\" type=\"Int32\">";
          for(int i = 0; i < NumberOfCells; i++){
            OutputFile << "  " << OutElementContainer[i].N1;
            OutputFile << "  " << OutElementContainer[i].N2;
            OutputFile << "  " << OutElementContainer[i].N3;
            OutputFile << "  " << OutElementContainer[i].N4;
          }
          OutputFile << "</DataArray>" << std::endl;

          OutputFile << "<DataArray Name=\"offsets\" format=\"ascii\" type=\"Int32\">";
          for(int i = 0; i < NumberOfCells; i++) OutputFile << "  " << (i+1)*4;
          OutputFile << "</DataArray>" << std::endl;


          OutputFile << "<DataArray Name=\"types\" format=\"ascii\" type=\"UInt8\">";
          for(int i = 0; i < NumberOfCells; i++) OutputFile << "  " << 9;
          OutputFile << "</DataArray>" << std::endl;

      OutputFile << "</Cells>" << std::endl;

      // Write Point Data
      OutputFile << "<PointData>" << std::endl;
        // Write Particle Volume
        OutputFile << "<DataArray Name=\"Mass\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
        for (auto &Node : OutNodeContainer) {
          OutputFile << "  " << Node.Mass;
        }
        OutputFile << "</DataArray>" << std::endl;

        // Write Particle Velocity
        OutputFile << "<DataArray Name=\"V\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">";
        for(int i = 0; i < NumberOfNodes; i++){
          OutputFile << "  " << CuttOff(OutNodeContainer[i].V[0]);
          OutputFile << "  " << CuttOff(OutNodeContainer[i].V[1]);
          OutputFile << "  " << CuttOff(OutNodeContainer[i].V[2]);
        }
        OutputFile << "</DataArray>" << std::endl;

//------ Write Nodal Momentum
        OutputFile << "<DataArray Name=\"Momentum\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">";
        for (auto &Node : OutNodeContainer) {
          OutputFile << "  " << CuttOff(Node.Momentum[0]);
          OutputFile << "  " << CuttOff(Node.Momentum[1]);
          OutputFile << "  " << CuttOff(Node.Momentum[2]);
        }
        OutputFile << "</DataArray>" << std::endl;

//------ Write Nodal Momentum
        OutputFile << "<DataArray Name=\"InternalForce\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">";
        for (auto &Node : OutNodeContainer) {
          OutputFile << "  " << CuttOff(Node.InternalForce[0]);
          OutputFile << "  " << CuttOff(Node.InternalForce[1]);
          OutputFile << "  " << CuttOff(Node.InternalForce[2]);
        }
        OutputFile << "</DataArray>" << std::endl;


      OutputFile << "</PointData>" << std::endl;
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
  }

void MPMOutputVTK::VTUParticleExport(std::string FileName, std::vector<MPMParticle> &Particles, std::vector<std::string> ParticleOutputDataNames){

    int NumberOfParticles  = Particles.size();

    //open a file stream for output
    std::ofstream OutputFile;
    OutputFile.open(FileName, std::ios::out);
    // Write headder
    OutputFile << "<?xml version=\"1.0\" ?>" << std::endl;
    OutputFile << "<VTKFile byte_order=\"LittleEndian\" type=\"UnstructuredGrid\" version=\"0.1\">" << std::endl;
    OutputFile << "<UnstructuredGrid>" << std::endl;
    // Write Piece Headder
    OutputFile << "<Piece NumberOfCells=\"" << NumberOfParticles << "\" NumberOfPoints=\"" << NumberOfParticles << "\">" << std::endl;
    // Write Points
    OutputFile << "<Points>" << std::endl;
    OutputFile << "<DataArray NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">";
    for(int i = 0; i < NumberOfParticles; i++){
      OutputFile << "  " << CuttOff(Particles[i].X[0]) ;
      OutputFile << "  " << CuttOff(Particles[i].X[1]) ;
      OutputFile << "  " << CuttOff(Particles[i].X[2]) ;
    }
    OutputFile << "</DataArray>" << std::endl;
    OutputFile << "</Points>" << std::endl;
    // Write Cells
    OutputFile << "<Cells>" << std::endl;

    OutputFile << "<DataArray Name=\"connectivity\" format=\"ascii\" type=\"Int32\">";
    for(int i = 0; i < NumberOfParticles; i++) OutputFile << "  " << i ;
    OutputFile << "</DataArray>" << std::endl;

    OutputFile << "<DataArray Name=\"offsets\" format=\"ascii\" type=\"Int32\">";
    for(int i = 0; i < NumberOfParticles; i++) OutputFile << "  " << (i+1)*1;
    OutputFile << "</DataArray>" << std::endl;

    OutputFile << "<DataArray Name=\"types\" format=\"ascii\" type=\"UInt8\">";
    for(int i = 0; i < NumberOfParticles; i++) OutputFile << "  " << 1;
    OutputFile << "</DataArray>" << std::endl;

    OutputFile << "</Cells>" << std::endl;

    // Write Point Data
    OutputFile << "<PointData>" << std::endl;
    for (std::string &PostTask : ParticleOutputDataNames){
      // define dummy output symbols
      int TaskKey;
      std::array<int   , 9> IntOut;
      std::array<double, 9> DoubleOut;

      // initial call to get post type
      TaskKey = Particles[0].GetPost(PostTask, DoubleOut, IntOut);
      // Post-Processing Interface
      // int: 0 -> key not defined !
      if (TaskKey==0) break;
      // int: 1 -> scalar integer
      if (TaskKey==1) OutputFile << "<DataArray Name=\""<< PostTask <<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">";
      // int: 2 -> scalar double
      if (TaskKey==2) OutputFile << "<DataArray Name=\""<< PostTask <<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">";
      // int: 3 -> vector[3] double
      if (TaskKey==3) OutputFile << "<DataArray Name=\""<< PostTask <<"\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">";
      // int: 4 -> tensor[9] double
      if (TaskKey==4) OutputFile << "<DataArray Name=\""<< PostTask <<"\" NumberOfComponents=\"9\" format=\"ascii\" type=\"Float64\">";

      // Then write data
      for (auto &Particle : Particles){
        TaskKey = Particle.GetPost(PostTask, DoubleOut, IntOut);
        if (TaskKey==1) OutputFile << "  " << IntOut[0];
        if (TaskKey==2) OutputFile << "  " << CuttOff(DoubleOut[0]);
        if (TaskKey==3) { OutputFile << "  " << CuttOff(DoubleOut[0]); OutputFile << "  " << CuttOff(DoubleOut[1]); OutputFile << "  " << CuttOff(DoubleOut[2]);}
        if (TaskKey==4) { for (int i=0;i<9;i++){OutputFile << "  " << CuttOff(DoubleOut[i]);}}
      }
      // close data field
      OutputFile << "</DataArray>" << std::endl;
    }
    OutputFile << "</PointData>" << std::endl;

    // Write Cell Data
    OutputFile << "<CellData/>" << std::endl;
    // Write Piece foot
    OutputFile << "</Piece>" << std::endl;
    // Write foot
    OutputFile << "</UnstructuredGrid>" << std::endl;
    OutputFile << "</VTKFile>" << std::endl;
    //close the file stram
    OutputFile.close();
};

void MPMOutputVTK::SetOutput(std::string FileName, std::vector<MPMParticle> &OutParticleContainer, std::vector<std::string> OutputStrings){
  // if FileName is valid:
  ParticleOutputFileNames.push_back(FileName);

  // if particle pointer valid ...
  // dummy
  ParticleOutputParticleContainers.push_back(&OutParticleContainer);

  // if datanames valid ...
  ParticleOutputDataNames.push_back(OutputStrings);

  // set write counter for the just described file to zero
  ParticleOutputFileWriteCounter.push_back(0);

  // if call is successfull
  NoParticleOutputs++;

}

void MPMOutputVTK::WriteOutput(double &t){
  // Write Particle Output Files
  for (int i=0;i<NoParticleOutputs;i++){
    std::string OutputFile = ParticleOutputFileNames[i]+ "_" + std::to_string(ParticleOutputFileWriteCounter[i]) + ".vtu";

    //write the output file
    VTUParticleExport(OutputFile, *(ParticleOutputParticleContainers[i]), ParticleOutputDataNames[i]);

    // For each written outputfile add an entry for pvd
    // conversion t to time need care as t can be very small
    std::stringstream stream ;
    stream << t ;
    std::string TIME = stream.str() ;
    std::vector<std::string> NewFile = {OutputFile,TIME};
    OutputFileContainer.push_back(NewFile);

    // Increment counter
    ParticleOutputFileWriteCounter[i]++;

  }


}

void MPMOutputVTK::WritePVD(){
  // Write PVD file mapping all exported .vtu onto a time stamp
  std::ofstream PVDFile;
  PVDFile.open(SimulationName+".pvd", std::ios::out);
  // Write headder
  PVDFile << "<VTKFile type=\"Collection\">" << std::endl;
  PVDFile << "<Collection>" << std::endl;

  for (std::vector<std::string> &WrittenFile : OutputFileContainer){

    PVDFile << "<DataSet ";
    PVDFile << "timestep=\""<< WrittenFile[1] << "\"";
    PVDFile << " ";
    PVDFile << "file=\""<< WrittenFile[0] << "\"";
    PVDFile << "/>" << std::endl;

  }

  // Finish file
  PVDFile << "</Collection>" << std::endl;
  PVDFile << "</VTKFile>" << std::endl;
  PVDFile.close();
}
