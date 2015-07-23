#pragma once

// Switch of some warnings for the Pluto headers
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wvla"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Weffc++"
#include "PParticle.h"
#pragma GCC diagnostic pop

#include <ostream>
#include <vector>

std::ostream& operator <<(std::ostream& stream, const PParticle* p);

std::ostream& operator <<(std::ostream& stream, const PParticle& p);

/**
 * @brief Get the ParticleTypeName of a const Pluto Particle
 * PParticle::Name() can't be called on const objects. Use this instead.
 * @param p
 * @return name of the particle type (e.g. "eta")
 */
std::string GetParticleName(const PParticle* p);
std::string GetParticleName(const int id);

void PrintParticleTable(std::ostream& stream, const std::vector<const PParticle*> plist);
