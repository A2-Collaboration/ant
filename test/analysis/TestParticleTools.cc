#include "catch.hpp"

#include "base/ParticleType.h"
#include "analysis/utils/ParticleTools.h"
#include "analysis/utils/RootAddons.h"

#include <cassert>
#include <iostream>


using namespace std;
using namespace ant;
using namespace ant::analysis::utils;


TEST_CASE("ParticleTools: GetPlutoProduction", "[analysis]") {

    using Ch_t = ParticleTypeTreeDatabase::Channel;

    auto get = [] (Ch_t ch) {
        return ParticleTools::GetPlutoProduction(ParticleTypeTreeDatabase::Get(ch));
    };
    CHECK(get(Ch_t::TwoPi0_4g) == "pi0 [ g g ] pi0 [ g g ] p");
    CHECK(get(Ch_t::ThreePi0_6g) == "pi0 [ g g ] pi0 [ g g ] pi0 [ g g ] p");
    CHECK(get(Ch_t::TwoPi0_2ggEpEm) == "pi0 [ g g ] pi0 [ e- g e+ ] p");
    CHECK(get(Ch_t::ThreePi0_4ggEpEm) == "pi0 [ g g ] pi0 [ g g ] pi0 [ e- g e+ ] p");
    CHECK(get(Ch_t::Pi0Eta_4g) == "eta [ g g ] pi0 [ g g ] p");
    CHECK(get(Ch_t::Pi0Eta_Pi03Pi0_8g) == "eta [ pi0 [ g g ] pi0 [ g g ] pi0 [ g g ] ] pi0 [ g g ] p");
    CHECK(get(Ch_t::Pi0Eta_gEpEm2g) == "eta [ g g ] pi0 [ e- g e+ ] p");
    CHECK(get(Ch_t::Pi0PiPi_2gPiPi) == "pi0 [ g g ] pi- pi+ p");
    CHECK(get(Ch_t::TwoPi0PiPi_4gPiPi) == "pi0 [ g g ] pi0 [ g g ] pi- pi+ p");
    CHECK(get(Ch_t::SigmaPlusK0s_6g) == "K0S [ pi0 [ g g ] pi0 [ g g ] ] Sigma+ [ pi0 [ g g ] p ]");
    CHECK(get(Ch_t::Omega_gEta_3g) == "w [ eta [ g g ] g ] p");
    CHECK(get(Ch_t::EtaPrime_gOmega_ggPi0_4g) == "eta' [ w [ g pi0 [ g g ] ] g ] p");
    CHECK(get(Ch_t::gp_pPi0) == "pi0 p");
}



