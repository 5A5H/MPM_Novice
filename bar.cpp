// 2D Material Point Method

#include <MPM_Process.hpp>
#include <MPM_OutputVTK.hpp>
#include <MPM_Particle.hpp>
#include <MPM_GridNode.hpp>
#include <MPM_GridNodeBC.hpp>
#include <MPM_TimeTracker.hpp>
#include <MPM_GridElement.hpp>
#include <MPM_SHPQ4.hpp>
#include <MPM_Read.hpp>
#include <MPM_Material.hpp>
#include <MPM_AceMaterials.hpp>
#include <MPM_HF.hpp>

#include <math.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>



//declare function
bool PointInQ4(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]);

//----------------------------------- Global Variables ----------------------------------------------------------------
std::vector<MPMParticle> Particle;
std::vector<MPMGridNode> GridNode;
std::vector<MPMGridElement> GridElement;
std::vector<MPMMaterial> Material;
double MassTolerance = 10e-6;

//----------------------------------- Global Functions ----------------------------------------------------------------


void ParticlesToGrid(){
// Calculate the grid nodal mass and momentum by mapping the particle mass and momemntum to the corresponding grid nodes.

//Loop over all Particles
for (auto &Pt : Particle){
  // Current Particles Position
  double Xp[3]; Xp[0] = Pt.X[0]; Xp[1] = Pt.X[1]; Xp[2] = Pt.X[2];
  if (false) std::cout << "X: " << Xp[0] << ", " << Xp[1] << ", " << Xp[2] <<std::endl;

  // Current Particles Mass
  double Mp = Pt.Mass;
  if (false) std::cout << "Mass: " << Mp <<std::endl;

  // Current Particles Density
  double rp = Pt.Density;
  if (false) std::cout << "Density: " << rp <<std::endl;

  // Current Particles Body Force
  double bp[3]; bp[0] = Pt.b[0]; bp[1] = Pt.b[1]; bp[2] = Pt.b[2];
  if (false) std::cout << "Body Force: " << bp[0] << ", " << bp[1] << ", " << bp[2] <<std::endl;

  // Current Particles Velocity
  double Vp[3]; Vp[0] = Pt.V[0]; Vp[1] = Pt.V[1]; Vp[2] = Pt.V[2];
  if (false) std::cout << "V: " << Vp[0] << ", " << Vp[1] << ", " << Vp[2] <<std::endl;

  // Current Cauchy Stresses
  double Sigp[3][3];
  Sigp[0][0] = Pt.Sig[0][0];Sigp[0][1] = Pt.Sig[0][1];Sigp[0][2] = Pt.Sig[0][2];
  Sigp[1][0] = Pt.Sig[1][0];Sigp[1][1] = Pt.Sig[1][1];Sigp[1][2] = Pt.Sig[1][2];
  Sigp[2][0] = Pt.Sig[2][0];Sigp[2][1] = Pt.Sig[2][1];Sigp[2][2] = Pt.Sig[2][2];
  if (false) {
    std::cout << "Sig: " << Sigp[0][0] << ", " << Sigp[0][1] << ", " << Sigp[0][2] <<std::endl;
    std::cout << "     " << Sigp[1][0] << ", " << Sigp[1][1] << ", " << Sigp[1][2] <<std::endl;
    std::cout << "     " << Sigp[2][0] << ", " << Sigp[2][1] << ", " << Sigp[2][2] <<std::endl;
  }

  //Loop over all Elements
  for (auto &Elmt : GridElement){

    // Current Element Nodes
    int ni[4]; ni[0] = Elmt.N1; ni[1] = Elmt.N2; ni[2] = Elmt.N3; ni[3] = Elmt.N4;
    if (false) std::cout << "Nodes : " << ni[0] << ", " << ni[1] << ", " << ni[2] << ", " << ni[3] <<std::endl;

    // Current Element Nodas Coordinates
    double XI[4][3];
    for (int i=0;i<4;i++){
      for (int j=0;j<3;j++)
      XI[i][j] = GridNode[ni[i]].X[j];
    }
    double X1[3]={XI[0][0],XI[0][1],XI[0][2]};
    double X2[3]={XI[1][0],XI[1][1],XI[1][2]};
    double X3[3]={XI[2][0],XI[2][1],XI[2][2]};
    double X4[3]={XI[3][0],XI[3][1],XI[3][2]};
    if (false) {
      std::cout << "Nodes 1 X: " << X1[0] << ", " << X1[1] << ", " << X1[2] <<std::endl;
      std::cout << "Nodes 2 X: " << X2[0] << ", " << X2[1] << ", " << X2[2] <<std::endl;
      std::cout << "Nodes 3 X: " << X3[0] << ", " << X3[1] << ", " << X3[2] <<std::endl;
      std::cout << "Nodes 4 X: " << X4[0] << ", " << X4[1] << ", " << X4[2] <<std::endl;
    }

    // Reset Current associated element
    Pt.Elmt = -1;

    //Check if Particle at Pt is inside Element Elmt
    if (PointInQ4(X1,X2,X3,X4,Xp)){
      // Particls is inside Element
      if (false) std::cout << "Particle Inside Element!" << std::endl;

      //Set this element as current element for the particle
      Pt.Elmt = Elmt.ID;

      //Loop over Element Nodes
      for (int i=0;i<4;i++){

        //Evaluate Shape Function
        MPMSHPQ4 Q4SHP;
        Q4SHP.evaluate(X1,X2,X3,X4,Xp);
        double NIP = Q4SHP.SHP(i);
        double DNIP[3];
        for (int j=0;j<3;j++){
          DNIP[j] = Q4SHP.SHP(i,j);
        }

        //Calculate Nodal Mass
        GridNode[ni[i]].Mass += Mp * NIP;

        //Calculate Nodal Momentum
        for (int dim=0;dim<3;dim++){
          GridNode[ni[i]].Momentum[dim] += Mp * Vp[dim] * NIP;
        }

        //Calculate Internal Nodal Force
        for (int dim=0;dim<3;dim++){
          for (int j=0;j<3;j++){
              GridNode[ni[i]].Force[dim] -= (Mp/rp) * Sigp[dim][j] * DNIP[j];
          }
        }

        //Calculate External Nodal Force
        for (int dim=0;dim<3;dim++){
        GridNode[ni[i]].Force[dim] += Mp * NIP * bp[dim];
        }


      // End Node Loop
      }

      break;
    }

  // End Element Loop
  }

// End Particle Loop
}

};
void GridTimeIntegration(double &dt){
  for (auto &Node : GridNode) {
    // Integrate the grid nodal momentum
    Node.Momentum[0] += Node.InternalForce[0] * dt;
    Node.Momentum[1] += Node.InternalForce[1] * dt;
    Node.Momentum[2] += Node.InternalForce[2] * dt;
    // Calculate the grid nodal velocity
    if (Node.Mass > MassTolerance){
      Node.V[0] = Node.Momentum[0]/Node.Mass;
      Node.V[1] = Node.Momentum[1]/Node.Mass;
      Node.V[2] = Node.Momentum[2]/Node.Mass;
    } else {
      Node.V[0] = 0e0;
      Node.V[1] = 0e0;
      Node.V[2] = 0e0;
    }
  }
};
void GridToParticle(double &dt, MPMMaterial &Mate){
  //Loop over all Particles
  for (auto &Pt : Particle){
    // Current Particles Position
    double Xp[3]; Xp[0] = Pt.X[0]; Xp[1] = Pt.X[1]; Xp[2] = Pt.X[2];
    if (false) std::cout << "X: " << Xp[0] << ", " << Xp[1] << ", " << Xp[2] <<std::endl;

    // Reset Particles Velocity Gradient
    for (int i=0;i<3;i++){
      for (int j=0;j<3;j++){
        Pt.L[i][j] = 0e0;
      }
    }

    // Load outdated Particle deformation gradient (in advnce of update)
    double Fn[3][3];
    for (int i=0;i<3;i++){
      for (int j=0;j<3;j++){
        Fn[i][j] = Pt.F[i][j];
      }
    }

    // Currently associated Element with this particle
    int AE = Pt.Elmt;

    // Current Element Nodes
    int n[4] = { GridElement[AE].N1, GridElement[AE].N2, GridElement[AE].N3, GridElement[AE].N4 };

    // Current Element Nodas Coordinates
    double XI[4][3];
    for (int i=0;i<4;i++){
      for (int j=0;j<3;j++)
      XI[i][j] = GridNode[n[i]].X[j];
    }
    double X1[3]={XI[0][0],XI[0][1],XI[0][2]};
    double X2[3]={XI[1][0],XI[1][1],XI[1][2]};
    double X3[3]={XI[2][0],XI[2][1],XI[2][2]};
    double X4[3]={XI[3][0],XI[3][1],XI[3][2]};

    // Evaluate shape functions
    MPMSHPQ4 Q4SHP;
    Q4SHP.evaluate(X1,X2,X3,X4,Xp);


    //Loop over all nodes of the associated element
    for (int I=0;I<4;I++){

      // Nodal Mass
      double MI = GridNode[n[I]].Mass;

      // Nodal Momentum
      double PI[3];
      for (int i=0;i<3;i++) PI[i] = GridNode[n[I]].Momentum[i];

      // Nodal Velocity
      double VI[3];
      for (int i=0;i<3;i++) VI[i] = GridNode[n[I]].V[i];

      // Nodal Force
      double FI[3];
      for (int i=0;i<3;i++) FI[i] = GridNode[n[I]].Force[i];

      //Nodal Shape function values
      double NIP = Q4SHP.SHP(I);
      double DNIP[3];
      for (int j=0;j<3;j++){
        DNIP[j] = Q4SHP.SHP(I,j);
      }

      // Check for nodal mass
      if (MI > MassTolerance){

        // Update Particle Position
        for (int dim=0;dim<3;dim++){
          Pt.X[dim] += (dt/MI) * NIP * PI[dim];
        }

        // Update Particle Velocity
        for (int dim=0;dim<3;dim++){
          Pt.V[dim] += (dt/MI) * NIP * FI[dim];
        }

      }// End check for nodal mass

      // Compute current velocity gradient
      for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
          Pt.L[i][j] += DNIP[j] * VI[i];
        }
      }

      // Identity Tensor
      double Iden[3][3];
      for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
          Iden[i][j] = 0e0;
        }
        Iden[i][i] = 1e0;
      }


      // Update deformation gradient
      for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
          Pt.F[i][j] = 0e0;
          for (int k=0;k<3;k++){
            Pt.F[i][j] += ( Iden[i][k] + Pt.L[i][k] * dt ) * Fn[k][j];
          }
        }
      }


      // Update Particles Stresses
      double dummyh[20];
      std::vector<double> dummyio;
      Mate.GetStresses(Pt.F, dummyh, Pt.Sig, dummyio);

    }// End nodal loop


  }//End particle loop
};




//------------------------------------------ MAIN ---------------------------------------------------------------------
int main()
{
    std::cout << "_____________________Welcome to MPM2D!____________________\n";
    // Genrate The Time Tracker
    MPMTimeTracker MPMTimings;
    MPMTimings.SetTime("Program Start");

    double t0 = 0.0;    double tmax = 6.3;
    double dt = 10e-5;  double rho  = 100;
    int step = 0;

//------------------------------------------- Output declaration ------------------------------------------------------
    bool ParaviewOutput = true;
    std::string ParticleOutputFile = "/Users/sash/mpm_2d/data/Bar/Bar_Particle";
    std::string GridOutputFile = "/Users/sash/mpm_2d/data/Bar/Bar_Grid";
    int PostFrequency = 200;
    MPMOutputVTK VTK1;
    MPM::File_CVS MyFile_V("/Users/sash/mpm_2d/data/Bar/MyCVS_V");
    MPM::File_CVS MyFile_X("/Users/sash/mpm_2d/data/Bar/MyCVS_X");

//------------------------------------------ Material declaration -----------------------------------------------------
    MPMMaterial Steel(3);
    double Emod = 10000;
    double nu   = 0.3;
    Steel.SetMaterialParameter(Emod);
    Steel.SetMaterialParameter(nu);
    Material.push_back(Steel);

//------------------------------------------ spatial discretization ---------------------------------------------------
    std::string InputfileParticle = "/Users/sash/mpm_2d/data/falling_bar_particle.txt";
    std::string InputfileNodes = "/Users/sash/mpm_2d/data/falling_bar_node.txt";
    std::string InputfileGrid = "/Users/sash/mpm_2d/data/falling_bar_element.txt";


    //Read and Create Objects
    {
    MPMTimings.SetTime("Read Start");
    ReadParticle(InputfileParticle, Particle);
    ReadGridNodes(InputfileNodes, GridNode);
    ReadGridElementsQ4(InputfileGrid, GridElement);
    std::cout << "Problem Data: " << std::endl;
    std::cout << "Number of Particles    : " << Particle.size() << std::endl;
    std::cout << "Number of Grid Nodes   : " << GridNode.size() << std::endl;
    std::cout << "Number of Grid Elements: " << GridElement.size() << std::endl;
    MPMTimings.SetTime("Read End");
    }

    //Set Particle Mass [and Partical initial condition experimental]
    MPMGridNodeBC MyFirstGridNodeBC;
    MyFirstGridNodeBC.setBC("EssentialBC","V",0);

    std::cout << "- Set Initial conditions" << std::endl;
    for (auto &Pt : Particle) {
       Pt.Mass = rho*Pt.Vol;
       Pt.V[1] = -0.1;
    }
     // Search and add nodes for bc
     for (auto &Node : GridNode) {
        if (Node.X[0]==0 || Node.X[0]==1 || Node.X[1]==0 || Node.X[1]==1) {
          MyFirstGridNodeBC.addGridNode(Node.ID);
        }
      }
//---------------------------------------------------------------------------------------------------------------------
// Check Materials
     for (auto &Mat : Material) {
        Mat.Report();
      }

//---------------------------------------------------------------------------------------------------------------------
    std::cout << "- Begin Time Integration" << std::endl;
    std::vector<int> PGC;  // PGC ->  ParticleGridConnectivity; holds element connectivity {0->element2 1->element3 .. noparticles->element89}
    std::string statusbar;
    MPMTimings.SetTime("Start TimeLoop");
    // Output initial VTK
    if (ParaviewOutput){
    VTK1.TestVTUGridExport(GridOutputFile + "_" + std::to_string(step) + ".vtu",GridNode,GridElement);
    VTK1.TestVTUParticleExport(ParticleOutputFile + "_" + std::to_string(step) + ".vtu",Particle);
    }
    // Time Loop
    if (false){
    for (double t=t0;t<tmax;t=t+dt){


    // Check for nan
    for (auto &Pt : Particle) {
        if (Pt.checkNAN()) return 1;
    }

    // Reset Grid
    for (auto &Node : GridNode) {
       Node.Reset();
    }

    // Find initial Element Connectivity
    //MPMTimings.SetTime("Search Start");
    //std::cout << "Start Initial Particle-Grid connectivity search   : " << std::endl;
    PGC.clear();
    for (auto &Pt : Particle) {
      PGC.push_back(-1);
      for (auto &Elmt : GridElement) {
         bool InsideThisElement;
         InsideThisElement = PointInQ4( GridNode[Elmt.N1].X, GridNode[Elmt.N2].X, GridNode[Elmt.N3].X, GridNode[Elmt.N4].X, Pt.X );
         if (InsideThisElement) {
           PGC.back() = Elmt.ID;
           break;
         }
       }
     }
     if (PGC.size() < Particle.size()){
       std::cout << "warning";
     }
    //MPMTimings.SetTime("Search End");

    // MPM Project Particle to Grid
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
        GridNode[PtElmt.N1].Momentum[i] += SHP.N1 * Pt.V[i] * Pt.Mass;
        GridNode[PtElmt.N2].Momentum[i] += SHP.N2 * Pt.V[i] * Pt.Mass;
        GridNode[PtElmt.N3].Momentum[i] += SHP.N3 * Pt.V[i] * Pt.Mass;
        GridNode[PtElmt.N4].Momentum[i] += SHP.N4 * Pt.V[i] * Pt.Mass;
        }
        // update nodal internal force vector
        double Sig11 = Pt.Stress[0];
        double Sig22 = Pt.Stress[1];
        double Sig12 = Pt.Stress[2];

        GridNode[PtElmt.N1].InternalForce[0] -= Pt.Vol * (Sig11*SHP.dN1dX + Sig12*SHP.dN1dY);
        GridNode[PtElmt.N1].InternalForce[1] -= Pt.Vol * (Sig12*SHP.dN1dX + Sig22*SHP.dN1dY);
        GridNode[PtElmt.N1].InternalForce[2] -= 0.0;

        GridNode[PtElmt.N2].InternalForce[0] -= Pt.Vol * (Sig11*SHP.dN2dX + Sig12*SHP.dN2dY);
        GridNode[PtElmt.N2].InternalForce[1] -= Pt.Vol * (Sig12*SHP.dN2dX + Sig22*SHP.dN2dY);
        GridNode[PtElmt.N2].InternalForce[2] -= 0.0;

        GridNode[PtElmt.N3].InternalForce[0] -= Pt.Vol * (Sig11*SHP.dN3dX + Sig12*SHP.dN3dY);
        GridNode[PtElmt.N3].InternalForce[1] -= Pt.Vol * (Sig12*SHP.dN3dX + Sig22*SHP.dN3dY);
        GridNode[PtElmt.N3].InternalForce[2] -= 0.0;

        GridNode[PtElmt.N4].InternalForce[0] -= Pt.Vol * (Sig11*SHP.dN4dX + Sig12*SHP.dN4dY);
        GridNode[PtElmt.N4].InternalForce[1] -= Pt.Vol * (Sig12*SHP.dN4dX + Sig22*SHP.dN4dY);
        GridNode[PtElmt.N4].InternalForce[2] -= 0.0;

      }
    }

    // Apply BC on grid nodal mom.
    // Time Integration
    for (auto &Node : GridNode) {
      if (Node.X[0]<10e-5 || Node.X[0]==1 || Node.X[1]<10e-5 || Node.X[1]==1) {
      //Node.Momentum[0] = 0;
      Node.Momentum[1] = 0;
      //Node.Momentum[2] = 0;
      //Node.InternalForce[0] = 0;
      Node.InternalForce[1] = 0;
      //Node.InternalForce[2] = 0;
    }
    }


    // Time Integration
    for (auto &Node : GridNode) {
      // time integrate momentum
      Node.Momentum[0] += Node.InternalForce[0] * dt;
      Node.Momentum[1] += Node.InternalForce[1] * dt;
      Node.Momentum[2] += Node.InternalForce[2] * dt;
      // time integrate velocity
      if (Node.Mass > MassTolerance){
        Node.V[0] = Node.Momentum[0]/Node.Mass;
        Node.V[1] = Node.Momentum[1]/Node.Mass;
        Node.V[2] = Node.Momentum[2]/Node.Mass;
      } else {
        Node.V[0] = 0.0;
        Node.V[1] = 0.0;
        Node.V[2] = 0.0;
      }
    }

    // Apply Boundary Conditions
    //MyFirstGridNodeBC.applyBC(GridNode);

    // MPM Project Grid to Particle
    for (auto &Pt : Particle) {
       //Element where the particle maps to
      if (PGC[Pt.ID]>=0){ //Catch node that is not in an element
        auto &PtElmt = GridElement[PGC[Pt.ID]];
        // Evaluate shape function
        SHP.evaluate(GridNode[PtElmt.N1].X, GridNode[PtElmt.N2].X, GridNode[PtElmt.N3].X, GridNode[PtElmt.N4].X, Pt.X);
        // update particle velocity and momentum

          if (GridNode[PtElmt.N1].Mass > MassTolerance) {
            for (int j=0;j<3;j++){
              Pt.V[j] += SHP.N1 * (GridNode[PtElmt.N1].InternalForce[j]/GridNode[PtElmt.N1].Mass)*dt;
              Pt.X[j] += SHP.N1 * (GridNode[PtElmt.N1].Momentum[j]/GridNode[PtElmt.N1].Mass)*dt;
            }
          }

          if (GridNode[PtElmt.N2].Mass > MassTolerance) {
            for (int j=0;j<3;j++){
              Pt.V[j] += SHP.N2 * (GridNode[PtElmt.N2].InternalForce[j]/GridNode[PtElmt.N2].Mass)*dt;
              Pt.X[j] += SHP.N2 * (GridNode[PtElmt.N2].Momentum[j]     /GridNode[PtElmt.N2].Mass)*dt;
            }
          }

          if (GridNode[PtElmt.N3].Mass > MassTolerance) {
            for (int j=0;j<3;j++){
              // !!!! EXAMPLE: double mass = .... and then use mass instead of ...
              Pt.V[j] += SHP.N3 * (GridNode[PtElmt.N3].InternalForce[j]/GridNode[PtElmt.N3].Mass)*dt;
              Pt.X[j] += SHP.N3 * (GridNode[PtElmt.N3].Momentum[j]     /GridNode[PtElmt.N3].Mass)*dt;
            }
          }

          if (GridNode[PtElmt.N4].Mass > MassTolerance) {
            for (int j=0;j<3;j++){
              Pt.V[j] += SHP.N4 * (GridNode[PtElmt.N4].InternalForce[j]/GridNode[PtElmt.N4].Mass)*dt;
              Pt.X[j] += SHP.N4 * (GridNode[PtElmt.N4].Momentum[j]     /GridNode[PtElmt.N4].Mass)*dt;
            }
          }

        // update particle deformation and stresses
        double Lp[4]; // 2D velocity gradient Lp = [ dvxdx , dvxdy ,dvydx, dvydy]
        Lp[0] = SHP.dN1dX * GridNode[PtElmt.N1].V[0] + SHP.dN2dX * GridNode[PtElmt.N2].V[0] + SHP.dN3dX * GridNode[PtElmt.N3].V[0] + SHP.dN4dX * GridNode[PtElmt.N4].V[0];
        Lp[1] = SHP.dN1dY * GridNode[PtElmt.N1].V[0] + SHP.dN2dY * GridNode[PtElmt.N2].V[0] + SHP.dN3dY * GridNode[PtElmt.N3].V[0] + SHP.dN4dY * GridNode[PtElmt.N4].V[0];
        Lp[2] = SHP.dN1dX * GridNode[PtElmt.N1].V[1] + SHP.dN2dX * GridNode[PtElmt.N2].V[1] + SHP.dN3dX * GridNode[PtElmt.N3].V[1] + SHP.dN4dX * GridNode[PtElmt.N4].V[1];
        Lp[3] = SHP.dN1dY * GridNode[PtElmt.N1].V[1] + SHP.dN2dY * GridNode[PtElmt.N2].V[1] + SHP.dN3dY * GridNode[PtElmt.N3].V[1] + SHP.dN4dY * GridNode[PtElmt.N4].V[1];
        double Fn[4] = { Pt.Deformation[0], Pt.Deformation[1], Pt.Deformation[2], Pt.Deformation[3]}; // 2D deformation gradient Fn = [ dxxdx , dxydx ,dxxdy, dxydy]
        //      update deformation gradient F = [ dxxdx , dxydx ,dxxdy, dxydy]
        double F[4];
        F[0] = Fn[0]*(1.0+Lp[0]*dt)+dt*Fn[2]*Lp[1];
        F[1] = Fn[1]*(1.0+Lp[0]*dt)+dt*Fn[3]*Lp[1];
        F[2] = Fn[2]*(1.0+Lp[3]*dt)+dt*Fn[0]*Lp[2];
        F[3] = Fn[3]*(1.0+Lp[3]*dt)+dt*Fn[1]*Lp[2];
        // Call Material Class
        double *SigMate;
        double defgrad[9]; //F = [F11,  F12,  F13,  F21,  F22,  F23,  F31,  F32,  F33]
        defgrad[0] = F[0];
        defgrad[1] = F[1];
        defgrad[3] = F[2];
        defgrad[4] = F[3];
        SigMate = Steel.getStresses(defgrad); //Sig = [Sig11,  Sig12,  Sig13,  Sig21,  Sig22,  Sig23,  Sig31,  Sig32,  Sig33]
        //      update stress and deformation onto particle
        Pt.Deformation[0] = F[0];
        Pt.Deformation[1] = F[1];
        Pt.Deformation[2] = F[2];
        Pt.Deformation[3] = F[3];
        Pt.Stress[0] = SigMate[0]; // Sig11
        Pt.Stress[1] = SigMate[4]; // Sig22
        Pt.Stress[2] = SigMate[1]; // Sig12

      }
    }

    // PostProcessing and Report
    // Export smth.
    *MyFile_V.strm() << t << " , " << Particle[0].V[0] << " , " << Particle[0].V[1] << " , " << Particle[0].V[2] << std::endl;
    *MyFile_X.strm() << t << " , " << Particle[0].X[0] << " , " << Particle[0].X[1] << " , " << Particle[0].X[2] << std::endl;
    // regular export
    if(step % PostFrequency == 0){
      // Progress Bar
      std::cout << "   Time Integration Progress : [";
      for (int i = 0;i<=(t/tmax)*24;i++) std::cout << "%";
      for (int i = 0;i<=24-(t/tmax)*24;i++) std::cout << "-";
      std::cout << "]";
      std::cout << "   Progress : " << std::setprecision(3) << std::setw(4) << std::left << (t/tmax)*100 << " % \r" << std::flush;
      // Paraview Output
      if (ParaviewOutput){
      VTK1.TestVTUGridExport(GridOutputFile + "_" + std::to_string(step) + ".vtu",GridNode,GridElement);
      VTK1.TestVTUParticleExport(ParticleOutputFile + "_" + std::to_string(step) + ".vtu",Particle);
      }
    }

    step++;
    }// end time loop
    }

    //New Loop

    // Check for nan
    for (auto &Pt : Particle) {
        if (Pt.checkNAN()) return 1;
    }

    // Reset Grid
    for (auto &Node : GridNode) {
       Node.Reset();
    }

    ParticlesToGrid();
    GridTimeIntegration(dt);
    GridToParticle(dt, Steel);


  std::cout << std::endl;
  MPMTimings.SetTime("End TimeLoop");
  std::cout << "- End Time Integration" << std::endl;


  MPMTimings.printTimeTable();
  std::cout << "_________________________ The End ________________________\n";
  MPMTimings.SetTime("Program Finish");
  return 0;
}

inline bool PointInQ4(double X1[3], double X2[3], double X3[3], double X4[3], double XP[3]){
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





// Some NOtES
// for (auto &Node : GlobalGridNodeContainer) {
//   std::cout << *(Node.X) << std::endl;
// }
// MPMTimings.printTimeTable();
// for (int i = 0; i < 10; i++) {
//         std::cout << "Status: " << i << "\r" << std::flush;
//         sleep(1);
// }
// std::cout << "Completed.\n";


// MPMTimings.SetTime("TestVTUExport Start");
// TestVTUGridExport("/Users/sash/mpm_2d/data/Grid_001.vtu",GridNode,GridElement);
// TestVTUParticleExport("/Users/sash/mpm_2d/data/Particle_001.vtu",Particle);
// MPMTimings.SetTime("TestVTUExport Finish");


// find  Particle
// double min=100;
//     int id;
//     for (auto &Pt : Particle) {
//        if (Pt.X[1] > min) {id = Pt.ID;}
//     }
//     std::cout << id <<std::endl;
//     Particle[id].Report();
