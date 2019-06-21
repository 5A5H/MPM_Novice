#include <LinearElasticity_A.hpp>
#include <string>
#include <map>
#include <array>

#include <ELSE_ContiMechLibrary.hpp>
#include <ELSE_System.hpp>

/*
Imolementation of a Linear Elasticity Law
Files:
      materials/LinearElasticity_A.hpp
      materials/LinearElasticity_A.cpp

  Following ELSE standard notation for tensor - vector transition:
    S[6] = [S_11, S_12, S_13, S_22, F_23, S_33]
    F[9] = [F_11, F_12, F_13, F_21, F_22, F_23, F_31, F_32, F_33]
*/

namespace ELSE{
namespace MPM{

// Constructor (baseically calls the base class Constructor)
LinearElasticity_A::LinearElasticity_A(std::string MaterialName):Material(MaterialName){};
// Destructor
LinearElasticity_A::~LinearElasticity_A(){};

// Actual Implementation of the constitutive law
int LinearElasticity_A::getCauchyStress(  std::array<double, 6> &CauchyStress,
                      std::array<double, 9> &DeformationGradient,
                      std::map<std::string, double> &MaterialHistory,
                      std::map<std::string, int>    &IntegerMaterialIO,
                      std::map<std::string, double> &DoubleMaterialIO){

                        // compute lame constants
                        double Emod,nue,lambda,mue;
                        getMaterialParameter("Emod",Emod);
                        getMaterialParameter("nue" ,nue );
                        ELSE::Conti::HookeToLame(Emod, nue, lambda, mue);

                        // compute linear stress tensor Eps_ij = 1/2 (F_ij + F_ji - 2 I_ij)
                        std::array<double, 6> Eps = {0,0,0,   0,0,     0};
                        std::array<double, 9> I   = {1,0,0, 0,1,0, 0,0,1};
                        // Eps_11, Eps_22, Eps_33
                        Eps[0] = 0.5e0 * (DeformationGradient[0] + DeformationGradient[0] - 2e0 * I[0]);
                        Eps[3] = 0.5e0 * (DeformationGradient[4] + DeformationGradient[4] - 2e0 * I[4]);
                        Eps[5] = 0.5e0 * (DeformationGradient[8] + DeformationGradient[8] - 2e0 * I[8]);
                        // Eps_12, Eps_13, Eps_23
                        Eps[1] = 0.5e0 * (DeformationGradient[1] + DeformationGradient[3]);
                        Eps[2] = 0.5e0 * (DeformationGradient[2] + DeformationGradient[6]);
                        Eps[4] = 0.5e0 * (DeformationGradient[5] + DeformationGradient[7]);

                        // compute trace of linear stresses
                        // trEps = Eps_11 + Eps_22 + Eps_33;
                        double trEps = Eps[0] + Eps[3] + Eps[5];

                        // compute linear stress tensor
                        // Sig_ij = lam * trEps * I_ij + 2*mue * Eps_ij
                        CauchyStress[0] = 2e0 * mue * Eps[0] + lambda * trEps;
                        CauchyStress[1] = 2e0 * mue * Eps[1];
                        CauchyStress[2] = 2e0 * mue * Eps[2];
                        CauchyStress[3] = 2e0 * mue * Eps[3] + lambda * trEps;
                        CauchyStress[4] = 2e0 * mue * Eps[4];
                        CauchyStress[5] = 2e0 * mue * Eps[5] + lambda * trEps;


                        //checkMaterialParameters(); // here only as example (should be checked only once)
                        return 0;
                      };

// Check wether all needed Material Parameters are present.
int LinearElasticity_A::checkMaterialParameters(){std::cout << "Checked Parametrs for LinearElasticity_A !" << std::endl; return 0;};


}
}
