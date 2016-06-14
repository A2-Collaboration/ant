#pragma once

#include "PParticle.h"

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
