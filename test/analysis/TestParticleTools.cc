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
    auto test = [] (std::function<std::string(const ParticleTypeTree&)> parser, const ParticleTypeTreeDatabase::Channel& channel, const std::string& result)
    {
        return parser(ParticleTypeTreeDatabase::Get(channel)) == result;
    };


    REQUIRE(test(ParticleTools::GetPlutoProduction,
                 ParticleTypeTreeDatabase::Channel::SigmaPlusK0s_6g,
                 "K0S [ pi0 [ g g ] pi0 [ g g ] ] Sigma+ [ pi0 [ g g ] p ] "));
    REQUIRE(test(ParticleTools::GetPlutoString,
                 ParticleTypeTreeDatabase::Channel::SigmaPlusK0s_6g,
                 "g p [ K0S [ pi0 [ g g ] pi0 [ g g ] ] Sigma+ [ pi0 [ g g ] p ] ] "));


    REQUIRE(test(ParticleTools::GetPlutoProduction,
                 ParticleTypeTreeDatabase::Channel::Omega_gEta_3g,
                 "omega [ eta [ g g ] g ] p "));
    REQUIRE(test(ParticleTools::GetPlutoString,
                 ParticleTypeTreeDatabase::Channel::Omega_gEta_3g,
                 "g p [ omega [ eta [ g g ] g ] p ] "));


    REQUIRE(test(ParticleTools::GetPlutoProduction,
                 ParticleTypeTreeDatabase::Channel::gp_pPi0,
                 "pi0 p "));

    REQUIRE(test(ParticleTools::GetPlutoString,
                 ParticleTypeTreeDatabase::Channel::gp_pPi0,
                 "g p [ pi0 p ] "));
}



