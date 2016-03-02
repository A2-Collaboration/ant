#include "catch.hpp"

#include "base/ParticleType.h"

using namespace ant;

TEST_CASE("ParticleType: Comparison", "[base]") {
    REQUIRE(ParticleTypeDatabase::BeamNeutron == ParticleTypeDatabase::BeamTarget);
    REQUIRE(ParticleTypeDatabase::BeamProton == ParticleTypeDatabase::BeamTarget);
    REQUIRE(ParticleTypeDatabase::BeamProton != ParticleTypeDatabase::BeamNeutron);

    REQUIRE(ParticleTypeDatabase::eCharged == ParticleTypeDatabase::eMinus);
    REQUIRE(ParticleTypeDatabase::eCharged == ParticleTypeDatabase::ePlus);
    REQUIRE(ParticleTypeDatabase::eMinus != ParticleTypeDatabase::ePlus);

    REQUIRE(ParticleTypeDatabase::PiCharged == ParticleTypeDatabase::PiMinus);
    REQUIRE(ParticleTypeDatabase::PiCharged == ParticleTypeDatabase::PiPlus);
    REQUIRE(ParticleTypeDatabase::PiMinus != ParticleTypeDatabase::PiPlus);
}


TEST_CASE("ParticleType: Photoproduction Threshold", "[base]") {
    REQUIRE(ParticleTypeDatabase::EtaPrime.PhotoproductionThresh() == Approx(1445.6).epsilon(0.1));
    REQUIRE(ParticleTypeDatabase::Pi0.PhotoproductionThresh() == Approx(144.7).epsilon(0.1));
}