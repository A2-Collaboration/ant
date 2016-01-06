#include "catch.hpp"

#include "base/ParticleType.h"
#include "analysis/data/Particle.h"
#include "analysis/utils/ParticleID.h"
#include "analysis/utils/root-addons.h"

#include "TCutG.h"

#include <cassert>
#include <iostream>


using namespace std;
using namespace ant;
using namespace ant::analysis::data;
using namespace ant::analysis::utils;

bool dEEtest(const CandidatePtr& cand, const std::shared_ptr<TCutG>& cut) {
    return cut->IsInside(cand->CaloEnergy, cand->VetoEnergy);
}

bool toftest(const CandidatePtr& cand, const std::shared_ptr<TCutG>& cut) {
    return cut->IsInside(cand->CaloEnergy, cand->Time);
}

void test_makeTCutG();

TEST_CASE("ParticleID:  makeTCutG", "[analysis]") {
    test_makeTCutG();
}


void test_nocuts();
void test_protoncut();
void test_electoncut();
void test_electonantprotoncut();
void test_tof();
void test_tofdee();


struct testdata {
    testdata();
    std::shared_ptr<TCutG> dEE_electron;
    std::shared_ptr<TCutG> dEE_proton;
    std::shared_ptr<TCutG> tofcut;

    CandidatePtr gamma;
    CandidatePtr neutron;
    CandidatePtr proton;
    CandidatePtr electron;
};

const testdata data;

TEST_CASE("ParticleID: NoCuts", "[analysis]") {
    test_nocuts();
}

TEST_CASE("ParticleID: proton cut", "[analysis]") {
    test_protoncut();
}

TEST_CASE("ParticleID: electron cut", "[analysis]") {
    test_electoncut();
}

TEST_CASE("ParticleID: electron + proton cut", "[analysis]") {
    test_electonantprotoncut();
}

TEST_CASE("ParticleID: tof cut", "[analysis]") {
    test_tof();
}

TEST_CASE("ParticleID: tof + dEE cut", "[analysis]") {
    test_tofdee();
}

void test_makeTCutG() {
    auto cut = root::makeTCutG("a", {{1,1},{3,1},{3,3},{1,3}});
    REQUIRE(cut->IsInside(2,2));
    REQUIRE_FALSE(cut->IsInside(2,4));

    auto cut2 = root::makeTCutG("b",{{1,1},{3,1},{3,3},{1,3}});
}

void test_nocuts() {
    BasicParticleID pid;

    // No cuts set
    REQUIRE(pid.Identify(data.gamma)   == &ParticleTypeDatabase::Photon);
    REQUIRE(pid.Identify(data.proton)  == nullptr); // no way to tell
    REQUIRE(pid.Identify(data.neutron) == &ParticleTypeDatabase::Photon); //No ToF/Size cut -> Neutrons will be photons
    REQUIRE(pid.Identify(data.electron)== nullptr); // no way to tell
}

void test_protoncut() {
    BasicParticleID pid;

    // proton cut
    pid.dEE_proton = data.dEE_proton;
    REQUIRE(pid.Identify(data.gamma)   == &ParticleTypeDatabase::Photon);
    REQUIRE(pid.Identify(data.proton)  == &ParticleTypeDatabase::Proton);
    REQUIRE(pid.Identify(data.neutron) == &ParticleTypeDatabase::Photon);
    REQUIRE(pid.Identify(data.electron)== nullptr);
}

void test_electoncut() {
    BasicParticleID pid;

    // electron cut
    pid.dEE_electron = data.dEE_electron;
    REQUIRE(pid.Identify(data.gamma)   == &ParticleTypeDatabase::Photon);
    REQUIRE(pid.Identify(data.proton)  == nullptr);
    REQUIRE(pid.Identify(data.neutron) == &ParticleTypeDatabase::Photon);
    REQUIRE(pid.Identify(data.electron)== &ParticleTypeDatabase::eCharged);

}

void test_electonantprotoncut() {
    BasicParticleID pid;

    // electron and proton cut
    pid.dEE_electron = data.dEE_electron;
    pid.dEE_proton = data.dEE_proton;
    REQUIRE(pid.Identify(data.gamma)   == &ParticleTypeDatabase::Photon);
    REQUIRE(pid.Identify(data.proton)  == &ParticleTypeDatabase::Proton);
    REQUIRE(pid.Identify(data.neutron) == &ParticleTypeDatabase::Photon);
    REQUIRE(pid.Identify(data.electron)== &ParticleTypeDatabase::eCharged);

}

void test_tof() {
    BasicParticleID pid;

    // hadronic/tof
    pid.tof = data.tofcut;
    REQUIRE(pid.Identify(data.gamma)   == &ParticleTypeDatabase::Photon);
    REQUIRE(pid.Identify(data.proton)  == &ParticleTypeDatabase::Proton);
    REQUIRE(pid.Identify(data.neutron) == &ParticleTypeDatabase::Neutron);
}


void test_tofdee() {
    BasicParticleID pid;

    // tof + dEE
    pid.dEE_electron = data.dEE_electron;
    pid.dEE_proton = data.dEE_proton;
    pid.tof = data.tofcut;
    REQUIRE(pid.Identify(data.gamma)   == &ParticleTypeDatabase::Photon);
    REQUIRE(pid.Identify(data.proton)  == &ParticleTypeDatabase::Proton);
    REQUIRE(pid.Identify(data.neutron) == &ParticleTypeDatabase::Neutron);
    REQUIRE(pid.Identify(data.electron)== &ParticleTypeDatabase::eCharged);

}




testdata::testdata()
{
    // set up cuts
    // Cut shapes used here are not realistic. Just some polygons for testing
    dEE_electron = root::makeTCutG("electron", {{50,2},{300,2},{300,2.5},{50,2.5}});
    dEE_proton = root::makeTCutG("proton",{{50,4},{300,4},{51,16}});
    tofcut = root::makeTCutG("tof",{{0,5},{400,5},{400,15},{0,15}});

    auto make_candidate = [] (double caloE, double time, double vetoE) {
        return std::make_shared<TCandidate>(Detector_t::Type_t::CB,
                                            caloE, 0, 0, time, 0, vetoE, 0, std::vector<TCluster>{});
    };

    //create particles and make sure they fulfill their respective cuts
    gamma     = make_candidate(100, 0, 0);
    assert(!dEEtest(gamma, dEE_proton));
    assert(!dEEtest(gamma, dEE_electron));
    assert(!toftest(gamma, tofcut));

    neutron   = make_candidate(100, 10, 0);
    assert(toftest(neutron,tofcut));
    assert(!dEEtest(neutron,dEE_proton));
    assert(!dEEtest(neutron,dEE_electron));

    proton    = make_candidate(100, 10, 10);
    assert(toftest(proton,tofcut));
    assert(dEEtest(proton,dEE_proton));
    assert(!dEEtest(proton,dEE_electron));

    electron  = make_candidate(100, 2, 2.1);
    assert(!toftest(electron,tofcut));
    assert(!dEEtest(electron,dEE_proton));
    assert(dEEtest(electron,dEE_electron));
}
