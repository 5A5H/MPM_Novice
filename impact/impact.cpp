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

#include <ELSE_GeometricLibrary.hpp>
#include <ELSE_ContiMechLibrary.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>


//declare function

void SetInitialCondition(std::vector<MPMParticle> &ParticleContainer, double rho, double V[3]){
  for (auto &Pt : ParticleContainer) {
     Pt.Mass = rho*Pt.Vol;
     Pt.Density = rho;
     Pt.V[0] = V[0]; Pt.V[1] = V[1]; Pt.V[2] = V[2];
   }
};
bool NextStep(std::vector<MPMParticle> &Tool, std::vector<MPMParticle> &Piece, std::vector<MPMGridNode> &GridNode, std::vector<MPMGridElement> &GridElement, double t){
  // Check for nan Tool
  for (auto &Pt : Tool) {
      if (Pt.checkNAN()) {
        std::cout << "NaN at Time :" << t << std::endl;
        Pt.Report();
        return true;
      }
  }
  // Check for nan Piece
  for (auto &Pt : Piece) {
      if (Pt.checkNAN()) {
        std::cout << "NaN at Time :" << t << std::endl;
        Pt.Report();
        return true;
      }
  }

  // Reset Grid
  for (auto &Node : GridNode) {
     Node.Reset();
  }
return false;
};
void ParticlesToGrid(std::vector<MPMParticle> &Particle, std::vector<MPMGridNode> &GridNode, std::vector<MPMGridElement> &GridElement){
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
    if (ELSE::Geometric::PointInQ4(X1,X2,X3,X4,Xp)){
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
void GridTimeIntegration(double &dt, std::vector<MPMGridNode> &GridNode, double &MassTolerance){
  for (auto &Node : GridNode) {
    // Integrate the grid nodal momentum
    Node.Momentum[0] += Node.Force[0] * dt;
    Node.Momentum[1] += Node.Force[1] * dt;
    Node.Momentum[2] += Node.Force[2] * dt;
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
    //Node.Report();
  }
};
void GridBoundaryCondition(std::vector<MPMGridNode> &GridNode){
  for (auto &Node : GridNode) {
     if (Node.X[0]<0.001) {
       Node.Momentum[0] = 0e0;
       Node.Force[0] = 0e0;
       Node.Momentum[1] = 0e0;
       Node.Force[1] = 0e0;
     }
     if (Node.X[1]==0) {
       Node.Momentum[1] = 0e0;
       Node.Force[1] = 0e0;
       Node.Momentum[0] = 0e0;
       Node.Force[0] = 0e0;
     }
     if (Node.X[0]>0.060-0.001) {
       Node.Momentum[0] = 0e0;
       Node.Force[0] = 0e0;
       Node.Momentum[1] = 0e0;
       Node.Force[1] = 0e0;
     }
   }
};
void GridToParticle(std::vector<MPMParticle> &Particle, std::vector<MPMGridNode> &GridNode, std::vector<MPMGridElement> &GridElement, double &dt, MPMMaterial &Mate, double &MassTolerance){
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

    // Check if node is in element
    if (AE >= 0) {

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
      if (false) std::cout <<"Mass :" << MI << std::endl;

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

    }// End nodal loop

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
      Mate.GetStresses(Pt.F, Pt.h, Pt.Sig, Pt.MateData);


  } // end if particle is alone in the dark ...
  }//End particle loop
};
void MoveRigidBody(double dt, std::vector<MPMParticle> &Particle){
  for (auto &Pt : Particle){
    Pt.V[0] = -.1;
    Pt.X[0] += -.1*dt;
  }
}


//----------------------------------- Global Variables ----------------------------------------------------------------
static std::vector<MPMParticle> Impactor,Target;
static std::vector<MPMGridNode> GridNode;
static std::vector<MPMGridElement> GridElement;
static double MassTolerance = 10e-6;

double t0 = 0.0; double tmax = 40.0e-6; double dt = 1.18e-8; int step = 0;
//double t0 = 0.0; double tmax = 40.0e-6; double dt = 1.18e-6; int step = 0;

bool ParaviewOutput = true; int PostFrequency = 10;
//std::string ImpactorOutputFile  = "/Users/sash/mpm_2d/impact/post/Impactor";
//std::string TargetOutputFile = "/Users/sash/mpm_2d/impact/post/Target";
//std::string GridOutputFile  = "/Users/sash/mpm_2d/impact/post/Grid";

// Setup Material Impactor
static int    MatIDs  = 5;
static double Emods   = 200.0e9;
static double nus     = 0.3;

// Setup Material Target
static int    MatID   = 6;
static double Emod    = 78.2e9;
static double nu      = 0.3;
static double y_0     = 300.0e6;
static double y_inf   = 10e10;
static double kh      = 0;
static double deltah  = 0;


std::string TargetInputFile       = "/Users/sash/mpm_2d/impact/Target.cvs";
std::string ImpactorInputFile     = "/Users/sash/mpm_2d/impact/Impactor.cvs";
std::string GridNodesInputFile    = "/Users/sash/mpm_2d/impact/Nodes.cvs";
std::string GridElementInputFile  = "/Users/sash/mpm_2d/impact/Elements.cvs";





int main(){
  std::cout << "_____________________Welcome to MPM2D!____________________\n";

  MPMOutputVTK VTKExport;
  MPMOutputVTK VTKOut("Impact");
  VTKOut.SetOutput("/Users/sash/mpm_2d/impact/post/TwoDisks_Plastic_Target", Target, {"SigMises","MaterialState","V"});
  VTKOut.SetOutput("/Users/sash/mpm_2d/impact/post/TwoDisks_Plastic_Impactor", Impactor, {});
  VTKOut.SetOutput("/Users/sash/mpm_2d/impact/post/TwoDisks_Plastic_Grid", GridNode, GridElement, {});


  MPMMaterial Steel(MatIDs);
  Steel.SetMaterialParameter(Emods);
  Steel.SetMaterialParameter(nus);

  MPMMaterial Alu(MatID);
  Alu.SetMaterialParameter(Emod);
  Alu.SetMaterialParameter(nu);
  Alu.SetMaterialParameter(y_0);
  Alu.SetMaterialParameter(y_inf);
  Alu.SetMaterialParameter(kh);
  Alu.SetMaterialParameter(deltah);

  ReadParticle(ImpactorInputFile, Impactor);
  ReadParticle(TargetInputFile, Target);
  ReadGridNodes(GridNodesInputFile, GridNode, 1);
  ReadGridElementsQ4(GridElementInputFile, GridElement);
  std::cout << "Problem Data: " << std::endl;
  std::cout << "Number of Particles    : " << Impactor.size() << std::endl;
  std::cout << "Number of Particles    : " << Target.size() << std::endl;
  std::cout << "Number of Grid Nodes   : " << GridNode.size() << std::endl;
  std::cout << "Number of Grid Elements: " << GridElement.size() << std::endl;

  double V0[3] = {0,-1160.0,0};
  SetInitialCondition(Impactor, 7850.0*2, V0);
  V0[0]=0e0; V0[1]=0e0; V0[2]=0e0;
  SetInitialCondition(Target, 2700.0, V0);
  // initial displace 0.03
  for (auto &Pt : Impactor) {Pt.X[1] -= 0.003;}

  VTKOut.WriteOutput(t0);


  //VTKExport.TestVTUParticleExport(ImpactorOutputFile  + "_" + std::to_string(step) + ".vtu", Impactor);
  //VTKExport.TestVTUParticleExport(TargetOutputFile    + "_" + std::to_string(step) + ".vtu", Target);
  //VTKExport.TestVTUGridExport(    GridOutputFile      + "_" + std::to_string(step) + ".vtu",GridNode,GridElement);

  for (double t=t0;t<tmax;t=t+dt){

    if (NextStep(Target, Impactor, GridNode, GridElement, t)) return 1;

    ParticlesToGrid(Target, GridNode, GridElement);
    ParticlesToGrid(Impactor, GridNode, GridElement);

    GridBoundaryCondition(GridNode);
    GridTimeIntegration(dt, GridNode, MassTolerance);

    GridToParticle(Target, GridNode, GridElement, dt, Alu, MassTolerance);
    GridToParticle(Impactor, GridNode, GridElement, dt, Steel, MassTolerance);

    // double maxCstress = 0.0;
    // double maxKstress = 0.0;
    // for (auto &Pt : Target){
    //   double sigmises;
    //   ELSE::Conti::VonMisesStress(Pt.Sig,sigmises);
    //   if (sigmises > 3e8){
    //     std::cout << "----------" << std::endl;
    //     std::cout << "Cauchy stress     : " << sigmises << std::endl;
    //     std::cout << "Kirchhoff stress  : " << Pt.MateData[3] << std::endl;
    //     std::cout << "MaterialState     : " << Pt.MateData[0] << std::endl;
    //     std::cout << "MaterialStatus    : " << Pt.MateData[1] << std::endl;
    //     std::cout << "MaterialIterations: " << Pt.MateData[2] << std::endl;
    //     //std::cout << "Jacobi            : " << Pt.MateData[2] << std::endl;
    //   }
    // }

    MPM::StatusBar(step, t, tmax, dt);
    if ((step % PostFrequency)==0) {
      VTKOut.WriteOutput(t);
      // double maxCstress = 0.0;
      // double maxKstress = 0.0;
      // for (auto &Pt : Target){
      //   double sigmises;
      //   ELSE::Conti::VonMisesStress(Pt.Sig,sigmises);
      //   maxCstress = max(maxCstress,sigmises);
      //   maxKstress = max(maxKstress,Pt.MateData[3]);
      // }
      // std::cout << "Max Cauchy stress: " << maxCstress << std::endl;
      // std::cout << "Max Kirchhoff stress: " << maxKstress << std::endl;
      // write out max SigMises and Taumises
      //VTKExport.TestVTUParticleExport(ImpactorOutputFile  + "_" + std::to_string(step) + ".vtu", Impactor);
      //VTKExport.TestVTUParticleExport(TargetOutputFile    + "_" + std::to_string(step) + ".vtu", Target);
      //VTKExport.TestVTUGridExport(    GridOutputFile      + "_" + std::to_string(step) + ".vtu",GridNode,GridElement);
    }
    step++;
  }
  std::cout << std::endl;
  VTKOut.WriteOutput(tmax);
  //VTKExport.TestVTUParticleExport(ImpactorOutputFile  + "_" + std::to_string(step) + ".vtu", Impactor);
  //VTKExport.TestVTUParticleExport(TargetOutputFile    + "_" + std::to_string(step) + ".vtu", Target);
  //VTKExport.TestVTUGridExport(    GridOutputFile      + "_" + std::to_string(step) + ".vtu",GridNode,GridElement);

  std::cout << "_________________________ The End ________________________\n";
  return 0;
}
