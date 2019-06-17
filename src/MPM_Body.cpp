#include <MPM_Body.hpp>

//Constructor
MPM_Body::MPM_Body(std::vector<MPMParticle> &ParticleContainer){
  // assign particle container with a pointer
  GlobalPartilceContainer=&ParticleContainer;
}

//Destructor
MPM_Body::~MPM_Body(){}
