#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "base/WrapTFile.h"
#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "analysis/input/pluto/PlutoReader.h"

#include "analysis/utils/Fitter.h"

#include "analysis/utils/MCFakeReconstructed.h"

#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;

void dotest_ideal(bool, bool);

TEST_CASE("Fitter: Ideal KinFitter, z vertex fixed, proton measured", "[analysis]") {
    dotest_ideal(false, false);
}

TEST_CASE("Fitter: Ideal KinFitter, z vertex fixed, proton UNmeasured", "[analysis]") {
    dotest_ideal(false, true);
}

TEST_CASE("Fitter: Ideal KinFitter, z vertex free, proton measured", "[analysis]") {
    dotest_ideal(true, false);
}

TEST_CASE("Fitter: Ideal KinFitter, z vertex free, proton UNmeasured", "[analysis]") {
    dotest_ideal(true, true);
}

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

void dotest_ideal(bool z_vertex, bool proton_unmeas) {

    auto rootfile = make_shared<WrapTFileInput>(string(TEST_BLOBS_DIRECTORY)+"/Pluto_Etap2g.root");
    PlutoReader reader(rootfile);

    REQUIRE_FALSE(reader.IsSource());

    utils::KinFitter kinfitter("kinfitter", 2,
                               make_shared<TestUncertaintyModel>(proton_unmeas),
                               z_vertex);
    if(z_vertex) {
        kinfitter.SetZVertexSigma(0);
        test::EnsureSetup(); // needed for MCFake
    }

    // use mc_fake with complete 4pi (no lost photons)
    auto mc_fake = z_vertex ? std_ext::make_unique<utils::MCFakeReconstructed>(true) : nullptr;

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
        auto& protons = eventdata.Particles.Get(ParticleTypeDatabase::Proton);
        auto& photons = eventdata.Particles.Get(ParticleTypeDatabase::Photon);

        REQUIRE(beam->Type() == ParticleTypeDatabase::BeamProton);
        REQUIRE(protons.size() == 1);
        REQUIRE(photons.size() == 2);

        kinfitter.SetEgammaBeam(beam->Ek());
        kinfitter.SetProton(protons.front());
        kinfitter.SetPhotons(photons);

        const APLCON::Result_t& res = kinfitter.DoFit();

        if(res.Status != APLCON::Result_Status_t::Success)
            continue;
        nFitOk++;

        REQUIRE(res.NIterations == 2);

        nFitIterations += res.NIterations;

        if(z_vertex)
            REQUIRE(kinfitter.GetFittedZVertex() == Approx(0.0).epsilon(1e-2)); // epsilon is quite high
    }

    CHECK(nEvents==500);
    CHECK(nFitOk==nEvents);
}
