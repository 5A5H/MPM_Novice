// 2D Material Point Method

#include <MPMProcess.hpp>
#include <MPMOutputVTK.hpp>
#include <MPM_Particle.hpp>
#include <MPM_GridNode.hpp>
#include <MPM_TimeTracker.hpp>
#include <MPM_GridElement.hpp>
#include <MPM_SHPQ4.hpp>

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

//declare function
void ReadParticleAreaData(std::vector<MPMParticle> &ParticleContainer);
void ReadParticlePosition(std::vector<MPMParticle> &ParticleContainer);
void ReadGridNodePosition(std::vector<MPMGridNode> &NodeContainer);
void ReadGridElementConnectivity(std::vector<MPMGridElement> &ElementContainer);
void TestVTUParticleExport(std::vector<MPMParticle> &OutParticleContainer);
void TestVTUGridExport(
  std::vector<MPMGridNode> &OutNodeContainer,
  std::vector<MPMGridElement> &OutElementContainer
);
bool PointInQ4(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);


//------------------------------------------ MAIN ---------------------------------------------
int main()
{
    std::cout << "_______Welcome to MPM2D!______\n";
    // Genrate The Time Tracker
    MPMTimeTracker MPMTimings;
    MPMTimings.SetTime("Program Start");

    // Base Data Containers for MPM
    std::vector<MPMParticle> GlobalParticleContainer;
    std::vector<MPMGridNode> GlobalGridNodeContainer;
    std::vector<MPMGridElement> GlobalGridElementContainer;

    // Some known Constants
    int NoGridNodes     = 576;
    int NoParticles     = 216;
    int NoGridElements  = 529;

    double Emod = 1000;
    double rho  = 1000;
    double nu   = 0.3;
    double lam  = (Emod*nu)/((1+nu)*(1-2*nu));
    double mue  = Emod/(2*(1+nu));

    //Create Objects
    for(int i=0; i<NoParticles;     i++) GlobalParticleContainer.push_back(   MPMParticle()   );
    for(int i=0; i<NoGridNodes;     i++) GlobalGridNodeContainer.push_back(   MPMGridNode()   );
    for(int i=0; i<NoGridElements;  i++) GlobalGridElementContainer.push_back(MPMGridElement());

    //Fill Objects
    ReadParticleAreaData(GlobalParticleContainer);
    ReadParticlePosition(GlobalParticleContainer);
    ReadGridNodePosition(GlobalGridNodeContainer);
    ReadGridElementConnectivity(GlobalGridElementContainer);

    //Set Particle Mass
    for(int i=0; i<NoParticles;     i++){
      GlobalParticleContainer[i].Density = rho;
      GlobalParticleContainer[i].Mass    = rho*GlobalParticleContainer[i].Vol;
    }

    // Find initial Element Connectivity
    // match elmt 48 part 1
    MPMTimings.SetTime("Search Start");
    int ParticleGridConnectivity[NoParticles];
    bool test;
    bool DetailedOutput = false;
    for(int i=0; i<GlobalParticleContainer.size(); i++){
      int TestElement;
      for(TestElement=0; TestElement<GlobalGridElementContainer.size(); TestElement++){
        if (PointInQ4(
          GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N1].X,
          GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N2].X,
          GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N3].X,
          GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N4].X,
          GlobalParticleContainer[i].X) ) {
          break;
        }
      }
    ParticleGridConnectivity[i] = TestElement;
    if (DetailedOutput) {
    std::cout << GlobalGridElementContainer[TestElement].ID << std::endl;
    std::cout << "Node1:" << GlobalGridElementContainer[TestElement].N1 << ": " << GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N1].X[0] << ", "<< GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N1].X[1] << ", "<< GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N1].X[2] << ", " << std::endl;
    std::cout << "Node2:" << GlobalGridElementContainer[TestElement].N2 << ": " << GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N2].X[0] << ", "<< GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N2].X[1] << ", "<< GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N2].X[2] << ", " << std::endl;
    std::cout << "Node3:" << GlobalGridElementContainer[TestElement].N3 << ": " << GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N3].X[0] << ", "<< GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N3].X[1] << ", "<< GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N3].X[2] << ", " << std::endl;
    std::cout << "Node4:" << GlobalGridElementContainer[TestElement].N4 << ": " << GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N4].X[0] << ", "<< GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N4].X[1] << ", "<< GlobalGridNodeContainer[GlobalGridElementContainer[TestElement].N4].X[2] << ", " << std::endl;
    }
    }
    MPMTimings.SetTime("Search End");

    // MPMSHPQ4 MyShape(
    //   GlobalGridNodeContainer[GlobalGridElementContainer[0].N1].X,
    //   GlobalGridNodeContainer[GlobalGridElementContainer[0].N2].X,
    //   GlobalGridNodeContainer[GlobalGridElementContainer[0].N3].X,
    //   GlobalGridNodeContainer[GlobalGridElementContainer[0].N4].X,
    //   GlobalParticleContainer[1].X
    // );

    // MPM Project Particle to grid

    MPMSHPQ4 SHP;
    for(int i=0; i<GlobalParticleContainer.size(); i++){
      // element for particle computation: ParticleGridConnectivity[i]
      // Cut Of if element is out of scope
      int ParticleElement = ParticleGridConnectivity[i];

      if (true) {

        // evaluate shape function
        SHP.evaluate(
          GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N1].X,
          GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N2].X,
          GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N3].X,
          GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N4].X,
          GlobalParticleContainer[i].X
        );

        // update nodal mass
        GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N1].Mass += SHP.N1 * GlobalParticleContainer[i].Mass;
        GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N2].Mass += SHP.N2 * GlobalParticleContainer[i].Mass;
        GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N3].Mass += SHP.N3 * GlobalParticleContainer[i].Mass;
        GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N4].Mass += SHP.N4 * GlobalParticleContainer[i].Mass;

      }// end if for cutoff criterion
    }// end loop over partcle


    for (auto &Node : GlobalGridNodeContainer) {
      std::cout << Node.Mass << std::endl;
    }


    //   GlobalParticleContainer.push_back(MPMParticle(i+1,Xp[i][0],Xp[i][1],0.0,Vp[i],dens)); // put new born particle in the global container
    //   GlobalParticleContainer[i].V[0] = (( Xp[i][0] < 0.5)? .1 : -0.1 );
    //   GlobalParticleContainer[i].V[1] = (( Xp[i][1] < 0.5)? .1 : -0.1 );
    //   GlobalParticleContainer[i].V[2] = 0.0;
    // }



    // First read in Particle Area data:
    // known particle len = 261
    //std::vector<double> Vp;
    //double dens = 10;
    // Read In The Problem:
              // Known number of particle


                 // Number of GridNodes

              // Number of GridElements
    //double Vp[NoParticles];         // Particle Volume Array
    //double Xp[NoParticles][2];      // Particle Position
    //double XI[NoGridNodes][2];         // Particle Position
    //int nElements[nElements][4];  // Grid Element Connectivity

    // ReadParticleAreaData(Vp);
    // ReadParticlePosition(Xp);
    // for(i=0;i<nParticle;i++){
    //   GlobalParticleContainer.push_back(MPMParticle(i+1,Xp[i][0],Xp[i][1],0.0,Vp[i],dens)); // put new born particle in the global container
    //   GlobalParticleContainer[i].V[0] = (( Xp[i][0] < 0.5)? .1 : -0.1 );
    //   GlobalParticleContainer[i].V[1] = (( Xp[i][1] < 0.5)? .1 : -0.1 );
    //   GlobalParticleContainer[i].V[2] = 0.0;
    // }
    // // Hover over all particles in the global container and get a report
    // std::vector<MPMParticle>::iterator v = GlobalParticleContainer.begin();
    // while( v != GlobalParticleContainer.end()) {
    // //std::cout << "value of volume = " << (*v).Volume << std::endl;
    // v++;
    // }
    // std::vector<MPMParticle>::iterator v = GlobalParticleContainer.begin();
    // while( v != GlobalParticleContainer.end()) {
    // (*v).Report();
    // v++;
    // }

    // double sumvol = 0;
    // std::vector<MPMParticle>::iterator v = GlobalParticleContainer.begin();
    // while( v != GlobalParticleContainer.end()) {
    // sumvol += (*v).Volume;
    // v++;
    // }
    // std::cout << "Complete Particle Volume :" << sumvol << std::endl;

    // code for look at an array
    //for(i=0;i<261;i++) std::cout << i << "  " << Vp[i] << std::endl;

    // MPMProcess MyProcess1; // Create a Process class (terminates automatically)
    // MPMOutputVTK MyOutput; // Create a MPMOutputVTK class (terminates automatically)
    // MyOutput.WriteVTK();

    // VTK Export
    MPMTimings.SetTime("TestVTUExport Start");
    TestVTUParticleExport(GlobalParticleContainer);
    TestVTUGridExport(GlobalGridNodeContainer,GlobalGridElementContainer);
    MPMTimings.SetTime("TestVTUExport Finish");


    std::cout << "___________ The End __________\n";
    MPMTimings.SetTime("Program Finish");
    MPMTimings.printTimeTable();
    return 0;
}

bool PointInQ4(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
  bool DetailedOutput = false;
  if (DetailedOutput){
  std::cout << "Test X1: " << X1[0] << " ," << X1[1] << " ," << X1[2] << std::endl;
  std::cout << "Test X2: " << X2[0] << " ," << X2[1] << " ," << X2[2] << std::endl;
  std::cout << "Test X3: " << X3[0] << " ," << X3[1] << " ," << X3[2] << std::endl;
  std::cout << "Test X4: " << X4[0] << " ," << X4[1] << " ," << X4[2] << std::endl;
  std::cout << "Test XP: " << XP[0] << " ," << XP[1] << " ," << XP[2] << std::endl;
  }

  // computation of tangents (in element order)
  double t[4][3];
  for(int j=0;j<3;j++){
      t[0][j] = X2[j]-X1[j];
      t[1][j] = X3[j]-X2[j];
      t[2][j] = X4[j]-X3[j];
      t[3][j] = X1[j]-X4[j];
  }
  if (DetailedOutput){
  std::cout << "TangentialVector X1X2: " << t[0][0] << ", " << t[0][1] << ", " << t[0][2] << ", " << std::endl;
  std::cout << "TangentialVector X2X3: " << t[1][0] << ", " << t[1][1] << ", " << t[1][2] << ", " << std::endl;
  std::cout << "TangentialVector X3X4: " << t[2][0] << ", " << t[2][1] << ", " << t[2][2] << ", " << std::endl;
  std::cout << "TangentialVector X4X1: " << t[3][0] << ", " << t[3][1] << ", " << t[3][2] << ", " << std::endl;
  }

  // compute normals
  double n[4][3];
  for(int j=0;j<4;j++){
  n[j][0] = -t[j][1]; n[j][1] = t[j][0]; n[j][2] = 0.0;
  }
  if (DetailedOutput){
  std::cout << "NormalVector X1X2: " << n[0][0] << ", " << n[0][1] << ", " << n[0][2] << ", " << std::endl;
  std::cout << "NormalVector X2X3: " << n[1][0] << ", " << n[1][1] << ", " << n[1][2] << ", " << std::endl;
  std::cout << "NormalVector X3X4: " << n[2][0] << ", " << n[2][1] << ", " << n[2][2] << ", " << std::endl;
  std::cout << "NormalVector X4X1: " << n[3][0] << ", " << n[3][1] << ", " << n[3][2] << ", " << std::endl;
  }
  // compute testvectors
  double tv[4][3];
  for(int j=0;j<3;j++){
      tv[0][j] = XP[j]-X1[j];
      tv[1][j] = XP[j]-X2[j];
      tv[2][j] = XP[j]-X3[j];
      tv[3][j] = XP[j]-X4[j];
  }
  if (DetailedOutput){
  std::cout << "TestVector X1XP: " << tv[0][0] << ", " << tv[0][1] << ", " << tv[0][2] << ", " << std::endl;
  std::cout << "TestVector X2XP: " << tv[1][0] << ", " << tv[1][1] << ", " << tv[1][2] << ", " << std::endl;
  std::cout << "TestVector X3XP: " << tv[2][0] << ", " << tv[2][1] << ", " << tv[2][2] << ", " << std::endl;
  std::cout << "TestVector X4XP: " << tv[3][0] << ", " << tv[3][1] << ", " << tv[3][2] << ", " << std::endl;
  }
  // compute sum of scalar products to comare with abs of scalar products
  double sumabs = 0.0;
  double sum    = 0.0;
  for(int i=0;i<4;i++){
  //std::cout << "ScalarProduct: " << ( tv[i][0]*n[i][0] + tv[i][1]*n[i][1] + tv[i][2]*n[i][2] ) << std::endl;
        sumabs  += abs( tv[i][0]*n[i][0] + tv[i][1]*n[i][1] + tv[i][2]*n[i][2] );
        sum     += ( tv[i][0]*n[i][0] + tv[i][1]*n[i][1] + tv[i][2]*n[i][2] );
  }
  if (DetailedOutput){
  std::cout << "sum : " << sum << std::endl;
  std::cout << "sumabs : " << sumabs << std::endl;
  }
  bool XPInside = (sumabs==sum)? true : false;
  return XPInside;
}

void ReadParticleAreaData(std::vector<MPMParticle> &ParticleContainer) {
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
      ParticleContainer[i].ID   = i;
      ParticleContainer[i].Vol  = std::stod(ParticleArea, &n);
      i++;
    } catch (const std::exception& e) {
      std::cout << "End Reading File" << std::endl;
    }
  }
  // Print for control the whole area
  //std::cout << "WholeArea " << WholeArea << std::endl;
}

void ReadParticlePosition(std::vector<MPMParticle> &ParticleContainer) {
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
      ParticleContainer[i].X[0] = std::stod(ParticleX, &nX);
      ParticleContainer[i].X[1] = std::stod(ParticleY, &nY);
      ParticleContainer[i].X[2] = 0.0;
      //vec[i][0] = std::stod(ParticleX, &nX);
      //vec[i][1] = std::stod(ParticleY, &nY);
      i++;
    } catch (const std::exception& e) {
      std::cout << "End Reading Particle Position File" << std::endl;
    }
  }
}

void ReadGridElementConnectivity(std::vector<MPMGridElement> &ElementContainer){
  // create file object
  std::ifstream ip("/Users/sash/mpm_2d/data/two_discs_element.txt"); // defiition of a ifstream object from the std space

  // check wether file could be loaded
  if(!ip.is_open()){
    std::cout << "ERROR: Cannot Open Particle Position File" << std::endl;
  } else {
    std::cout << "Start Reading Particle Position File" << std::endl;
  }

  //declare variables
  std::string n1;
  std::string n2;
  std::string n3;
  std::string n4;
  std::string::size_type nn1;
  std::string::size_type nn2;
  std::string::size_type nn3;
  std::string::size_type nn4;


  //loop over lines of file
  int i = 0;
  while(ip.good()){
    getline(ip, n1, ',');
    getline(ip, n2, ',');
    getline(ip, n3, ',');
    getline(ip, n4, '\n');
    try{
      ElementContainer[i].ID = i+1;
      ElementContainer[i].N1 = std::stoi(n1, &nn1);
      ElementContainer[i].N2 = std::stoi(n2, &nn2);
      ElementContainer[i].N3 = std::stoi(n3, &nn3);
      ElementContainer[i].N4 = std::stoi(n4, &nn4);
      //vec[i][0] = std::stod(ParticleX, &nX);
      //vec[i][1] = std::stod(ParticleY, &nY);
      i++;
    } catch (const std::exception& e) {
      std::cout << "End Reading Particle Position File" << std::endl;
    }
  }
}

void ReadGridNodePosition(std::vector<MPMGridNode> &NodeContainer){
  // create file object
  std::ifstream ip("/Users/sash/mpm_2d/data/two_discs_node.txt"); // defiition of a ifstream object from the std space

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
        NodeContainer[i].ID   = i;
        NodeContainer[i].X[0] = std::stod(ParticleX, &nX);
        NodeContainer[i].X[1] = std::stod(ParticleY, &nY);
        NodeContainer[i].X[2] = 0.0;
        //vec[i][0] = std::stod(ParticleX, &nX);
        //vec[i][1] = std::stod(ParticleY, &nY);
        i++;
      } catch (const std::exception& e) {
        std::cout << "End Reading Particle Position File" << std::endl;
      }
    }
}

void TestVTUParticleExport(std::vector<MPMParticle> &OutParticleContainer){


  //collecting information
  //piece information
  int NumberOfParticles  =   OutParticleContainer.size();
  int NumberOfCells   =   0;


  //open a file stream for output
  std::ofstream OutputFile;
  OutputFile.open("/Users/sash/mpm_2d/data/vpuout_2.vtu", std::ios::out);
  // Write headder
  OutputFile << "<?xml version=\"1.0\" ?>" << std::endl;
  OutputFile << "<VTKFile byte_order=\"LittleEndian\" type=\"UnstructuredGrid\" version=\"0.1\">" << std::endl;
  OutputFile << "<UnstructuredGrid>" << std::endl;

  // Piece 1 -> Particle
  // Write Piece Headder
  OutputFile << "<Piece NumberOfCells=\"" << NumberOfParticles << "\" NumberOfPoints=\"" << NumberOfParticles << "\">" << std::endl;
  // Write Points
  OutputFile << "<Points>" << std::endl;

  OutputFile << "<DataArray NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">";
  for(int i = 0; i < NumberOfParticles; i++){
    OutputFile << "  " << OutParticleContainer[i].X[0] ;
    OutputFile << "  " << OutParticleContainer[i].X[1] ;
    OutputFile << "  " << OutParticleContainer[i].X[2] ;
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
    OutputFile << "<DataArray Name=\"Volume\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\">";
    for(int i = 0; i < NumberOfParticles; i++){
      OutputFile << "  " << OutParticleContainer[i].Vol ;
    }
    OutputFile << "</DataArray>" << std::endl;

    // Write Particle Velocity
    OutputFile << "<DataArray Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">";
    for(int i = 0; i < NumberOfParticles; i++){
      OutputFile << "  " << OutParticleContainer[i].V[0] ;
      OutputFile << "  " << OutParticleContainer[i].V[1] ;
      OutputFile << "  " << OutParticleContainer[i].V[2] ;
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

void TestVTUGridExport(
  std::vector<MPMGridNode> &OutNodeContainer,
  std::vector<MPMGridElement> &OutElementContainer){
    {


      //collecting information
      //piece information
      int NumberOfNodes   =   OutNodeContainer.size();
      int NumberOfCells   =   OutElementContainer.size();


      //open a file stream for output
      std::ofstream OutputFile;
      OutputFile.open("/Users/sash/mpm_2d/data/vpuoutgrid_2.vtu", std::ios::out);
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
            OutputFile << "  " << OutElementContainer[i].N1-1 ;
            OutputFile << "  " << OutElementContainer[i].N2-1 ;
            OutputFile << "  " << OutElementContainer[i].N3-1 ;
            OutputFile << "  " << OutElementContainer[i].N4-1 ;
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
        OutputFile << "<DataArray Name=\"Mass\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\">";
        for (auto &Node : OutNodeContainer) {
          OutputFile << "  " << Node.Mass;
        }
        OutputFile << "</DataArray>" << std::endl;

        // Write Particle Velocity
        OutputFile << "<DataArray Name=\"V\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">";
        for(int i = 0; i < NumberOfNodes; i++){
          OutputFile << "  " << OutNodeContainer[i].V[0] ;
          OutputFile << "  " << OutNodeContainer[i].V[1] ;
          OutputFile << "  " << OutNodeContainer[i].V[2] ;
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
// Some NOtES
// for (auto &Node : GlobalGridNodeContainer) {
//   std::cout << *(Node.X) << std::endl;
// }
