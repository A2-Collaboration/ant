#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "base/WrapTFile.h"
#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "analysis/input/pluto/PlutoReader.h"

#include "analysis/utils/Fitter.h"

#include "analysis/utils/MCFakeReconstructed.h"
#include "analysis/utils/MCSmear.h"

#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;

void dotest(bool, bool, bool);


TEST_CASE("Fitter: Ideal KinFitter, z vertex fixed, proton measured", "[analysis]") {
    dotest(false, false, false);
}

TEST_CASE("Fitter: Ideal KinFitter, z vertex fixed, proton UNmeasured", "[analysis]") {
    dotest(false, true, false);
}

TEST_CASE("Fitter: Ideal KinFitter, z vertex free, proton measured", "[analysis]") {
    dotest(true, false, false);
}

TEST_CASE("Fitter: Ideal KinFitter, z vertex free, proton UNmeasured", "[analysis]") {
    dotest(true, true, false);
}

TEST_CASE("Fitter: Smeared KinFitter, z vertex fixed, proton measured", "[analysis]") {
    dotest(false, false, true);
}

TEST_CASE("Fitter: Smeared KinFitter, z vertex fixed, proton UNmeasured", "[analysis]") {
    dotest(false, true, true);
}

//TEST_CASE("Fitter: Smeared KinFitter, z vertex free, proton measured", "[analysis]") {
//    dotest(true, false, true);
//}

//TEST_CASE("Fitter: Smeared KinFitter, z vertex free, proton UNmeasured", "[analysis]") {
//    dotest(true, true, true);
//}

struct TestUncertaintyModel : utils::UncertaintyModel {
    const bool ProtonUnmeasured;
    TestUncertaintyModel(bool protonUnmeasured) : ProtonUnmeasured(protonUnmeasured) {}
    virtual utils::Uncertainties_t GetSigmas(const TParticle& particle) const
    {
        utils::Uncertainties_t  u{1, 0.1, 0.1}; // any value should work...
        if(ProtonUnmeasured && particle.Type() == ParticleTypeDatabase::Proton)
            u.sigmaE = 0;
        return u;
    }
};

void dotest(bool z_vertex, bool proton_unmeas, bool smeared) {

    auto rootfile = make_shared<WrapTFileInput>(string(TEST_BLOBS_DIRECTORY)+"/Pluto_Etap2g.root");
    PlutoReader reader(rootfile);

    REQUIRE_FALSE(reader.IsSource());

    auto model = make_shared<TestUncertaintyModel>(proton_unmeas);

    utils::KinFitter kinfitter("kinfitter", 2,
                               model, z_vertex);
    if(z_vertex) {
        kinfitter.SetZVertexSigma(5.0);
        test::EnsureSetup(); // needed for MCFake
    }

    // use mc_fake with complete 4pi (no lost photons)
    auto mc_fake = z_vertex ? std_ext::make_unique<utils::MCFakeReconstructed>(true) : nullptr;
    auto mc_smear = smeared ? std_ext::make_unique<utils::MCSmear>(model) : nullptr;

    unsigned nEvents = 0;
    unsigned nFitOk = 0;
    unsigned nFitIterations = 0;
    while(true) {
        TEvent event;
        if(!reader.ReadNextEvent(event))
            break;
        nEvents++;

        INFO("nEvents="+to_string(nEvents));

        const TEventData& eventdata = z_vertex ? mc_fake->Get(event.MCTrue()) : event.MCTrue();

        const TParticlePtr& beam = event.MCTrue().ParticleTree->Get();
        TParticleList protons = eventdata.Particles.Get(ParticleTypeDatabase::Proton);
        TParticleList photons = eventdata.Particles.Get(ParticleTypeDatabase::Photon);

        REQUIRE(beam->Type() == ParticleTypeDatabase::BeamProton);
        REQUIRE(protons.size() == 1);
        REQUIRE(photons.size() == 2);

        TParticlePtr proton = protons.front();

        if(smeared) {
            proton = mc_smear->Smear(proton);
            for(auto& photon : photons)
                photon = mc_smear->Smear(photon);
        }

        kinfitter.SetEgammaBeam(beam->Ek());
        kinfitter.SetProton(proton);
        kinfitter.SetPhotons(photons);

        const APLCON::Result_t& res = kinfitter.DoFit();

        if(res.Status != APLCON::Result_Status_t::Success)
            continue;
        nFitOk++;

        if(!smeared) // in ideal conditions, the fitter should converge immediately
            REQUIRE(res.NIterations == 2);

        nFitIterations += res.NIterations;

        if(z_vertex)
            REQUIRE(kinfitter.GetFittedZVertex() == Approx(0.0).epsilon(1e-3)); // epsilon is quite high
    }

    CHECK(nEvents==500);
    CHECK(nFitOk==nEvents);
    if(smeared)
        CHECK(nFitIterations > nEvents*2);
}
