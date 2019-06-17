// MPMBody
#ifndef _MPM_BODY_HPP
#define _MPM_BODY_HPP

#include <vector>
#include <MPM_Particle.hpp>

class MPM_Body {
  public:
      MPM_Body(std::vector<MPMParticle> &ParticleContainer);
      ~MPM_Body();

  private:

    // Pointer to global container for the actual particle objects this body consists of
    std::vector<MPMParticle> *GlobalPartilceContainer;

    //Indexes of belonging particles
    std::vector<int> BodyParticles;

};

#endif
