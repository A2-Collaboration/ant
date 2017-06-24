#pragma once

#include "ParticleType.h"
#include "base/std_ext/math.h"

namespace ant {
namespace math {


/**
 * @brief Calculate the center of mass energy W = sqrt(s) from beam photon energy and fixed target particle type
 * @param Eg Photon beam energy
 * @param target particle type
 * @return W = sqrt(s)
 */
inline double W(const double Eg, const ParticleTypeDatabase::Type& target) {
    return sqrt(std_ext::sqr(Eg+target.Mass())-std_ext::sqr(Eg));
}

}
}
