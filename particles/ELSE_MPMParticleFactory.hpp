// Definition of the MPM Particle Factory
#ifndef _ELSE_MPM_PARTICLE_FACTORY_HPP_
#define _ELSE_MPM_PARTICLE_FACTORY_HPP_

#include <string>
#include <map>
#include <array>

// Inlcude Particle Base Class
#include <ELSE_MPMParticle.hpp>


/*
The Particler Factory
  - The Factory is called to create particle objects.
  - It returns a pointer to these objects.
*/

namespace ELSE{
namespace MPM{

  Particle* CreateParticle(std::string ParticleKEY, int &ID, std::array<double, 3> &Position){
    // ParticleKEY is a string identifier to various types of Particles
    // ID is an integer identifier for the new particle
    // Position is a vector size array denoting the initial position of the new particle
    if (ParticleKEY.size()==0) return new Particle(ID, Position);

    // if function evaluates to here there is no implementation for the requested particle key
    std::cerr << "Error: No Implementation for a Particle: " << ParticleKEY << " found." << std::endl;
    return nullptr;
  }



}
}


#endif
