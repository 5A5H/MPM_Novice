// 2D Material Point Method

#include <MPMProcess.hpp>
#include <MPMOutputVTK.hpp>
#include <MPM_Particle.hpp>
#include <MPM_GridNode.hpp>
#include <MPM_TimeTracker.hpp>
#include <MPM_GridElement.hpp>
#include <MPM_SHPQ4.hpp>
#include <MPM_Read.hpp>

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

//declare function
void TestVTUParticleExport(std::string FileName, std::vector<MPMParticle> &OutParticleContainer);
void TestVTUGridExport(
  std::string FileName,
  std::vector<MPMGridNode> &OutNodeContainer,
  std::vector<MPMGridElement> &OutElementContainer
);
bool PointInQ4(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);
//----------------------------------- Global Variables ----------------------------------------
std::vector<MPMParticle> Particle;
std::vector<MPMGridNode> GridNode;
std::vector<MPMGridElement> GridElement;

//------------------------------------------ MAIN ---------------------------------------------
int main()
{
    std::cout << "______________Welcome to MPM2D!_____________\n";
    // Genrate The Time Tracker
    MPMTimeTracker MPMTimings;
    MPMTimings.SetTime("Program Start");

    double Emod = 1000;
    double rho  = 1000;
    double nu   = 0.3;
    double lam  = (Emod*nu)/((1+nu)*(1-2*nu));
    double mue  = Emod/(2*(1+nu));
    double dt   = 0.001;


    //Read and Create Objects
    MPMTimings.SetTime("Read Start");
    ReadParticle("/Users/sash/mpm_2d/data/two_discs_particledata.txt", Particle);
    ReadGridNodes("/Users/sash/mpm_2d/data/two_discs_node.txt", GridNode);
    ReadGridElementsQ4("/Users/sash/mpm_2d/data/two_discs_element.txt", GridElement);
    std::cout << "Problem Data: " << std::endl;
    std::cout << "Number of Particles    : " << Particle.size() << std::endl;
    std::cout << "Number of Grid Nodes   : " << GridNode.size() << std::endl;
    std::cout << "Number of Grid Elements: " << GridElement.size() << std::endl;
    MPMTimings.SetTime("Read End");

    //Set Particle Mass [and BC experimental]
    for (auto &Pt : Particle) {
       Pt.Mass = rho*Pt.Vol;
       if (Pt.X[0]<0.5){
         Pt.V[0] = 0.1; Pt.V[1] = 0.1; Pt.V[2] = 0.0;
       } else {
         Pt.V[0] = -0.1; Pt.V[1] = -0.1; Pt.V[2] = 0.0;
       }
     }

    // Find initial Element Connectivity
    MPMTimings.SetTime("Search Start");
    std::vector<int> PGC;  // PGC ->  ParticleGridConnectivity; holds element connectivity {0->element2 1->element3 .. noparticles->element89}
    for (auto &Pt : Particle) {
      for (auto &Elmt : GridElement) {
         bool InsideThisElement;
         InsideThisElement = PointInQ4( GridNode[Elmt.N1].X, GridNode[Elmt.N2].X, GridNode[Elmt.N3].X, GridNode[Elmt.N4].X, Pt.X );
         if (InsideThisElement) {
           PGC.push_back(Elmt.ID);
           break;
         }
       }
     }
    MPMTimings.SetTime("Search End");

    // MPM Project Particle to grid
    MPMSHPQ4 SHP;
    for (auto &Pt : Particle) {
       //Element where the particle maps to
      if (PGC[Pt.ID]>=0){ //Catch node that is not in an element
        auto &PtElmt = GridElement[PGC[Pt.ID]];
        // Evaluate shape function
        SHP.evaluate(GridNode[PtElmt.N1].X, GridNode[PtElmt.N2].X, GridNode[PtElmt.N3].X, GridNode[PtElmt.N4].X, Pt.X);
        // update nodal masses
        GridNode[PtElmt.N1].Mass += SHP.N1 * Pt.Mass;
        GridNode[PtElmt.N2].Mass += SHP.N2 * Pt.Mass;
        GridNode[PtElmt.N3].Mass += SHP.N3 * Pt.Mass;
        GridNode[PtElmt.N4].Mass += SHP.N4 * Pt.Mass;
        // update nodal momentum
        for (int i=0;i<3;i++) {
        GridNode[PtElmt.N1].Momentum[i] += SHP.N1 * Pt.V[0] * Pt.Mass;
        GridNode[PtElmt.N2].Momentum[i] += SHP.N2 * Pt.V[0] * Pt.Mass;
        GridNode[PtElmt.N3].Momentum[i] += SHP.N3 * Pt.V[0] * Pt.Mass;
        GridNode[PtElmt.N4].Momentum[i] += SHP.N4 * Pt.V[0] * Pt.Mass;
        }
        // update nodal internal force vector
        for (int i=0;i<3;i++) {
        GridNode[PtElmt.N1].InternalForce[i] += 0.0;
        GridNode[PtElmt.N2].InternalForce[i] += 0.0;
        GridNode[PtElmt.N3].InternalForce[i] += 0.0;
        GridNode[PtElmt.N4].InternalForce[i] += 0.0;
        }
      }
    }


    // Time Integration
    for (auto &Node : GridNode) {
      // time integrate momentum
      Node.Momentum[0] += Node.InternalForce[0] * dt;
      Node.Momentum[1] += Node.InternalForce[1] * dt;
      Node.Momentum[2] += Node.InternalForce[2] * dt;
      // time integrate velocity
      if (Node.Mass > 10e-6){
        Node.V[0] += Node.Momentum[0] * dt;
        Node.V[1] += Node.Momentum[1] * dt;
        Node.V[2] += Node.Momentum[2] * dt;
      } else {
        Node.V[0] = Node.Momentum[0]/Node.Mass;
        Node.V[1] = Node.Momentum[1]/Node.Mass;
        Node.V[2] = Node.Momentum[2]/Node.Mass;
      }

    }

    // for(int i=0; i<GlobalParticleContainer.size(); i++){
    //   // element for particle computation: ParticleGridConnectivity[i]
    //   // Cut Of if element is out of scope
    //   int ParticleElement = ParticleGridConnectivity[i];
    //
    //   if (true) {
    //
    //     // evaluate shape function
    //     SHP.evaluate(
    //       GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N1].X,
    //       GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N2].X,
    //       GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N3].X,
    //       GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N4].X,
    //       GlobalParticleContainer[i].X
    //     );
    //
    //     // update nodal mass
    //     GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N1].Mass += SHP.N1 * GlobalParticleContainer[i].Mass;
    //     GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N2].Mass += SHP.N2 * GlobalParticleContainer[i].Mass;
    //     GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N3].Mass += SHP.N3 * GlobalParticleContainer[i].Mass;
    //     GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N4].Mass += SHP.N4 * GlobalParticleContainer[i].Mass;
    //     if(SHP.N1!=SHP.N1 || SHP.N2!=SHP.N2 || SHP.N3!=SHP.N3 || SHP.N4!=SHP.N4){
    //       std::cout << SHP.N1 << "  ";
    //       std::cout << SHP.N2 << "  ";
    //       std::cout << SHP.N3 << "  ";
    //       std::cout << SHP.N4 << "  ";
    //       std::cout <<std::endl;
    //       std::cout << "shp evaluated for";
    //       std::cout <<std::endl;
    //       std::cout << "X1 :";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N1].X[0] << " ,";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N1].X[1] << " ,";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N1].X[2] << " ,";
    //       std::cout << "| X2 :";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N2].X[0] << " ,";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N2].X[1] << " ,";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N2].X[2] << " ,";
    //       std::cout << "| X3 :";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N3].X[0] << " ,";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N3].X[1] << " ,";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N3].X[2] << " ,";
    //       std::cout << "| X4 :";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N4].X[0] << " ,";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N4].X[1] << " ,";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N4].X[2] << " ,";
    //       std::cout <<std::endl;
    //       std::cout << "| XP :";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N4].X[0] << " ,";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N4].X[1] << " ,";
    //       std::cout << GlobalGridNodeContainer[GlobalGridElementContainer[ParticleElement].N4].X[2] << " ,";
    //       std::cout <<std::endl;
    //       std::cout << "The element is";
    //       std::cout << ParticleElement;
    //       std::cout << "  with id  ";
    //       std::cout << GlobalGridElementContainer[ParticleElement].ID;
    //       std::cout <<std::endl;
    //       GlobalGridElementContainer[ParticleElement].Report();
    //     }
    //
    //   }// end if for cutoff criterion
    // }// end loop over partcle


    // for (auto &Node : GlobalGridNodeContainer) {
    //   std::cout << Node.ID << Node.Mass << std::endl;
    // }
    // for (int i=0;i<GlobalParticleContainer.size();i++) {
    //   std::cout << ParticleGridConnectivity[i] << std::endl;
    // }

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
    // GlobalGridElementContainer[0].Report();
    // GlobalGridElementContainer[528].Report();
    // GlobalGridElementContainer[529].Report();
    // GlobalGridElementContainer[530].Report();
    // GlobalGridElementContainer[531].Report();
    // std::cout << "size " << GlobalGridElementContainer.size();
    // VTK Export
    MPMTimings.SetTime("TestVTUExport Start");
    TestVTUGridExport("/Users/sash/mpm_2d/data/Grid_001.vtu",GridNode,GridElement);
    TestVTUParticleExport("/Users/sash/mpm_2d/data/Particle_001.vtu",Particle);
    MPMTimings.SetTime("TestVTUExport Finish");


    std::cout << "__________________ The End _________________\n";
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
  bool XPInside = (sumabs-sum < 10e-10)? true : false;
  if (DetailedOutput){
  std::cout << "sum      : " << sum << std::endl;
  std::cout << "sumabs   : " << sumabs << std::endl;
  if(XPInside) {
      std::cout << "XPInside : True" << std::endl;
    } else {
      std::cout << "XPInside : False" << std::endl;
    }
  }
  return XPInside;
}

void TestVTUParticleExport(std::string FileName, std::vector<MPMParticle> &OutParticleContainer){


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
    // Write Particle Mass
    OutputFile << "<DataArray Name=\"Mass\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float32\">";
    for (auto &Particle : OutParticleContainer) {
      OutputFile << "  " << Particle.Mass;
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

//------ Write Nodal Momentum
        OutputFile << "<DataArray Name=\"Momentum\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">";
        for (auto &Node : OutNodeContainer) {
          OutputFile << "  " << Node.Momentum[0];
          OutputFile << "  " << Node.Momentum[1];
          OutputFile << "  " << Node.Momentum[2];
        }
        OutputFile << "</DataArray>" << std::endl;

//------ Write Nodal Momentum
        OutputFile << "<DataArray Name=\"InternalForce\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float32\">";
        for (auto &Node : OutNodeContainer) {
          OutputFile << "  " << Node.InternalForce[0];
          OutputFile << "  " << Node.InternalForce[1];
          OutputFile << "  " << Node.InternalForce[2];
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
