#ifndef _ELSE_MPM_PARTICLE_HPP_
#define _ELSE_MPM_PARTICLE_HPP_

#include <math.h>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <ELSE_ContiMechLibrary.hpp>

#include <ELSE_MPMMaterial.hpp>

/*
The MPM Particle Base Class
- a request by the material is made to store allocate an amout of data for the material history and io as well as their keys for post processing
- conatins many standard data, where an ending _t denotes quantities from the previous time step

Scalar Quantities:
 * ID -> unique id for the particle object
 * Mass -> mass associated with the particle (computed from volume and ref. density from material)
 * Volume -> reference volume (may change over time by detF)

Vector Quantities:
 * X -> reference position
 * x -> current position
 * v -> current velocity
 * ForceInt -> current internal forces resulting from material stresses
 * ForceExt -> current external forces, either from natural BC or from body forces

SymmetricTensor Quantities:
 * CauchyStress -> Cauchy stress tensor [S11, S12, S13, S22, S23, S33]

Tensor Quantities:
 * F -> deformation gradient [F11, F12, F13, F21, F22, F23, F31, F32, F33]
 * L -> spacial velocity gradient [L11, L12, L13, L21, L22, L23, L31, L32, L33]

Material Data at particle:
 * h, ht -> history field (size is controlled by material model)
 * IntMIO -> integer material io data
 * DblMIO -> double material io data



*/

using std::array;
using std::vector;
using std::map;
using std::string;

namespace ELSE {
namespace MPM {

class Particle{
public:

  // constructor
  Particle();
  Particle(int ID, array<double, 3> Position);
  Particle(int ID, array<double, 3> Position, double Volume);
  Particle(int ID, array<double, 3> Position, double Volume, Material *particlesmaterial);
  // destructor
  ~Particle();

  // material interaction
  void setMaterial(Material *particlesmaterial);

  // functions to interact with particle data
  void get(string KEY, int &Value);
  void get(string KEY, double &Value);
  void get(string KEY, array<double, 3> &Value);
  void get(string KEY, array<double, 6> &Value);
  void get(string KEY, array<double, 9> &Value);

  void set(string KEY, int &Value);
  void set(string KEY, double &Value);
  void set(string KEY, array<double, 3> &Value);
  void set(string KEY, array<double, 6> &Value);
  void set(string KEY, array<double, 9> &Value);

  void update(string KEY, int &Value);
  void update(string KEY, double &Value);
  void update(string KEY, array<double, 3> &Value);
  void update(string KEY, array<double, 6> &Value);
  void update(string KEY, array<double, 9> &Value);

  // // standard functions for functionality
  // virtual void updateStresses();

private:


  // for conveniend access to particle data prepare 4 map classes <-> the constructor maps properly
  map<string, int> ParticleIntData;
  map<string, double> ParticleDblData;
  map<string, array<double, 3>> ParticleVectorData;
  map<string, array<double, 6>> ParticleSymTensorData;
  map<string, array<double, 9>> ParticleTensorData;

  // for material interaction
  Material *ParticlesMaterial;
  map<string, int> IntMIO;
  map<string, double> h,h_t,DblMIO;

};

}
}

#endif
