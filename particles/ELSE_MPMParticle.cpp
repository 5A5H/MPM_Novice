#include <ELSE_MPMParticle.hpp>

namespace ELSE {
namespace MPM {

// constructors of the particle base class
// 1. Id and Position
Particle::Particle(int ID, array<double, 3> Position){
  ParticleIntData["ID"] = ID;
  ParticleVectorData["X"] = Position;
}
// 2. ID, Position and Volume
Particle::Particle(int ID, array<double, 3> Position, double Volume){
  ParticleIntData["ID"] = ID;
  ParticleVectorData["X"] = Position;
  ParticleDblData["V0"] = Volume;
}
// 3. ID, Position, Volume and Material
Particle::Particle(int ID, array<double, 3> Position, double Volume, Material *particlesmaterial){
  ParticleIntData["ID"] = ID;
  ParticleVectorData["X"] = Position;
  ParticleDblData["V0"] = Volume;
  ParticlesMaterial = particlesmaterial;
}

// destructor of the particle base class
Particle::~Particle(){}

// associate material to particle
void Particle::setMaterial(Material *particlesmaterial){ParticlesMaterial = particlesmaterial;};

// getter functions
void Particle::get(string KEY, int              &Value){Value = ParticleIntData.find(KEY) -> second;}
void Particle::get(string KEY, double           &Value){Value = ParticleDblData.find(KEY) -> second;}
void Particle::get(string KEY, array<double, 3> &Value){Value = ParticleVectorData.find(KEY) -> second;}
void Particle::get(string KEY, array<double, 6> &Value){Value = ParticleSymTensorData.find(KEY) -> second;}
void Particle::get(string KEY, array<double, 9> &Value){Value = ParticleTensorData.find(KEY) -> second;}

// setter functions
void Particle::set(string KEY, int              &Value){ParticleIntData[KEY]      = Value;}
void Particle::set(string KEY, double           &Value){ParticleDblData[KEY]      = Value;}
void Particle::set(string KEY, array<double, 3> &Value){ParticleVectorData[KEY]   = Value;}
void Particle::set(string KEY, array<double, 6> &Value){ParticleSymTensorData[KEY]= Value;}
void Particle::set(string KEY, array<double, 9> &Value){ParticleTensorData[KEY]   = Value;}

// increment functions
void Particle::update(string KEY, int              &Value){ParticleIntData[KEY]=ParticleIntData.find(KEY) -> second + Value;}
void Particle::update(string KEY, double           &Value){ParticleDblData[KEY]=ParticleDblData.find(KEY) -> second + Value;}
void Particle::update(string KEY, array<double, 3> &Value){
  array<double, 3> PValue = ParticleVectorData.find(KEY) -> second;
  for (int i=0; i<3; i++) PValue[i] += Value[i];
  ParticleVectorData[KEY]   = PValue;
}
void Particle::update(string KEY, array<double, 6> &Value){
  array<double, 6> PValue = ParticleSymTensorData.find(KEY) -> second;
  for (int i=0; i<6; i++) PValue[i] += Value[i];
  ParticleSymTensorData[KEY]   = PValue;
}
void Particle::update(string KEY, array<double, 9> &Value){
  array<double, 9> PValue = ParticleTensorData.find(KEY) -> second;
  for (int i=0; i<9; i++) PValue[i] += Value[i];
  ParticleTensorData[KEY]   = PValue;
}

// interface function updateStresses (is virtual)
void Particle::updateStress(){
  std::array<double, 6> Sig;
  std::array<double, 9> F;
  get("F",F);

  ParticlesMaterial -> getCauchyStress( Sig, F, h, IntMIO, DblMIO);

  set("CauchyStress",Sig);

};

}
}
