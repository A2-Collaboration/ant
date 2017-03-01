#include "catch.hpp"

#include "base/ParticleType.h"
#include "analysis/utils/particle_tools.h"
#include "analysis/utils/root-addons.h"

#include <cassert>
#include <iostream>


using namespace std;
using namespace ant;
using namespace ant::analysis::utils;


TEST_CASE("ParticleTools: GetPlutoProduction", "[analysis]") {

    CHECK(ParticleTools::GetPlutoProduction(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::SigmaPlusK0s_6g)) ==
         "K0S [ pi0 [ g g ] pi0 [ g g ] ] Sigma+ [ pi0 [ g g ] p ]");
    CHECK(ParticleTools::GetPlutoProduction(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Omega_gEta_3g)) ==
          "w [ eta [ g g ] g ] p");
    CHECK(ParticleTools::GetPlutoProduction(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::gp_pPi0)) ==
          "pi0 p");
}



