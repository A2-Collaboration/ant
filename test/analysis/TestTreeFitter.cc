#include "catch.hpp"
#include "catch_config.h"
#include "expconfig_helpers.h"

#include "base/WrapTFile.h"
#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "analysis/input/pluto/PlutoReader.h"

#include "analysis/utils/fitter/TreeFitter.h"

#include "analysis/utils/MCFakeReconstructed.h"
#include "analysis/utils/MCSmear.h"
#include "analysis/utils/ParticleTools.h"

#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;

void dotest_Etap2g();
void dotest_EtapOmegaG_simple();
void dotest_EtapOmegaG_filter(bool);

TEST_CASE("TreeFitter: Etap2g: NoFilter", "[analysis]") {
    dotest_Etap2g();
}

TEST_CASE("TreeFitter: EtapOmegaG: NoFilter", "[analysis]") {
    dotest_EtapOmegaG_simple();
}

TEST_CASE("TreeFitter: EtapOmegaG: Filter, cut", "[analysis]") {
    dotest_EtapOmegaG_filter(false);
}

TEST_CASE("TreeFitter: EtapOmegaG: Filter, sort", "[analysis]") {
    dotest_EtapOmegaG_filter(true);
}

struct TestUncertaintyModel : utils::UncertaintyModel {

    const utils::A2SimpleGeometry geo;

    TestUncertaintyModel() {}
    virtual utils::Uncertainties_t GetSigmas(const TParticle& particle) const
    {
        utils::Uncertainties_t  u{
                    0.05*particle.Ek(),
                    std_ext::degree_to_radian(2.0),
                    std_ext::degree_to_radian(2.0),
                    15 // shower depth in cm
        };
        auto detector = geo.DetectorFromAngles(particle);
        if(detector== Detector_t::Any_t::None)
            detector = Detector_t::Type_t::CB;

        if(detector & Detector_t::Type_t::CB) {
            u.sigmaCB_R = 0.5;
        }
        else if(detector & Detector_t::Type_t::TAPS) {
            u.sigmaTAPS_Rxy = 8;
            u.sigmaTAPS_L = 0.5;
        }

        if(particle.Type() == ParticleTypeDatabase::Proton)
            u.sigmaEk = 0;
        return u;
    }
};

void dotest_Etap2g() {
    test::EnsureSetup();

    auto rootfile = make_shared<WrapTFileInput>(string(TEST_BLOBS_DIRECTORY)+"/Pluto_Etap2g.root");
    PlutoReader reader(rootfile);

    auto model = make_shared<TestUncertaintyModel>();

    utils::TreeFitter treefitter(
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2g),
                model, true);

    treefitter.SetZVertexSigma(3.0);

    // use mc_fake with complete 4pi (no lost photons)
    utils::MCFakeReconstructed mc_fake(true);

    unsigned nEvents = 0;
    unsigned nFailed = 0;

    while(true) {
        input::event_t event;
        if(!reader.ReadNextEvent(event))
            break;
        nEvents++;

        INFO("nEvents="+to_string(nEvents));

        auto mctrue_particles = mc_fake.Get(event.MCTrue());

        TParticlePtr beam = event.MCTrue().ParticleTree->Get();
        TParticleList protons = mctrue_particles.Get(ParticleTypeDatabase::Proton);
        TParticleList photons = mctrue_particles.Get(ParticleTypeDatabase::Photon);

        REQUIRE(beam->Type() == ParticleTypeDatabase::BeamProton);
        REQUIRE(protons.size() == 1);
        REQUIRE(photons.size() == 2);

        TParticlePtr proton = protons.front();
        LorentzVec photon_sum;
        for(auto& photon : photons)
            photon_sum += *photon;
        LorentzVec constraint_unsmeared = *beam - *proton - photon_sum;

        REQUIRE(constraint_unsmeared.E == Approx(0).epsilon(1e-3));
        REQUIRE(constraint_unsmeared.p.x == Approx(0).epsilon(1e-3));
        REQUIRE(constraint_unsmeared.p.y == Approx(0).epsilon(1e-3));
        REQUIRE(constraint_unsmeared.p.z == Approx(0).epsilon(1e-3));


        // do the fit
        treefitter.PrepareFits(beam->Ek(), proton, photons);
        APLCON::Result_t res;

        unsigned nPerms = 0;
        double prb = std_ext::NaN;
        unsigned bestPerm = 0;
        while(treefitter.NextFit(res)) {
            nPerms++;
            if(res.Status != APLCON::Result_Status_t::Success)
                continue;
            if(!std_ext::copy_if_greater(prb, res.Probability))
                continue;
            bestPerm = nPerms;
        }
        REQUIRE(nPerms == 1);
        if(prb != Approx(1.0)) {
            nFailed++;
            continue;
        }
        REQUIRE(bestPerm == 1);

    }

    REQUIRE(nFailed == 28);
    REQUIRE(nEvents == 1000);
}

void dotest_EtapOmegaG_simple() {
    test::EnsureSetup();

    auto rootfile = make_shared<WrapTFileInput>(string(TEST_BLOBS_DIRECTORY)+"/Pluto_EtapOmegaG.root");
    PlutoReader reader(rootfile);

    auto model = make_shared<TestUncertaintyModel>();

    utils::TreeFitter treefitter(
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_gOmega_ggPi0_4g),
                model, true);

    treefitter.SetZVertexSigma(3.0);


    // use mc_fake with complete 4pi (no lost photons)
    utils::MCFakeReconstructed mc_fake(true);

    unsigned nEvents = 0;
    unsigned nFailed = 0;

    while(true) {
        input::event_t event;
        if(!reader.ReadNextEvent(event))
            break;
        nEvents++;

        INFO("nEvents="+to_string(nEvents));

        auto mctrue_particles = mc_fake.Get(event.MCTrue());

        TParticlePtr beam = event.MCTrue().ParticleTree->Get();
        TParticleList protons = mctrue_particles.Get(ParticleTypeDatabase::Proton);
        TParticleList photons = mctrue_particles.Get(ParticleTypeDatabase::Photon);

        REQUIRE(beam->Type() == ParticleTypeDatabase::BeamProton);
        REQUIRE(protons.size() == 1);
        REQUIRE(photons.size() == 4);

        TParticlePtr proton = protons.front();
        LorentzVec photon_sum;
        for(auto& photon : photons)
            photon_sum += *photon;
        LorentzVec constraint_unsmeared = *beam - *proton - photon_sum;

        REQUIRE(constraint_unsmeared.E == Approx(0).epsilon(1e-3));
        REQUIRE(constraint_unsmeared.p.x == Approx(0).epsilon(1e-3));
        REQUIRE(constraint_unsmeared.p.y == Approx(0).epsilon(1e-3));
        REQUIRE(constraint_unsmeared.p.z == Approx(0).epsilon(1e-3));


        // do the fit
        treefitter.PrepareFits(beam->Ek(), proton, photons);
        APLCON::Result_t res;

        unsigned nPerms = 0;
        double prb = std_ext::NaN;
        unsigned bestPerm = 0;
        while(treefitter.NextFit(res)) {
            nPerms++;
            if(res.Status != APLCON::Result_Status_t::Success)
                continue;
            if(!std_ext::copy_if_greater(prb, res.Probability))
                continue;
            bestPerm = nPerms;
        }
        REQUIRE(nPerms == 12);
        if(prb != Approx(1.0)) {
            nFailed++;
            continue;
        }
        REQUIRE(bestPerm == 1);

    }

    REQUIRE(nFailed == 3);
    REQUIRE(nEvents == 100);
}


void dotest_EtapOmegaG_filter(bool sort) {
    test::EnsureSetup();

    auto rootfile = make_shared<WrapTFileInput>(string(TEST_BLOBS_DIRECTORY)+"/Pluto_EtapOmegaG.root");
    PlutoReader reader(rootfile);

    auto model = make_shared<TestUncertaintyModel>();

    utils::TreeFitter treefitter(
                ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_gOmega_ggPi0_4g),
                model, true);

    treefitter.SetZVertexSigma(3.0);

    auto fitted_Pi0 = treefitter.GetTreeNode(ParticleTypeDatabase::Pi0);
    REQUIRE(fitted_Pi0);

    if(sort) {
        // sort by inverse chi2
        treefitter.SetIterationFilter([fitted_Pi0] () {
            auto& node = fitted_Pi0->Get();
            return 1.0/std_ext::sqr(ParticleTypeDatabase::Pi0.Mass() - node.LVSum.M());
        },
        3 // do the best 3 iterations...
        );
    }
    else {
        // simple cut
        treefitter.SetIterationFilter([fitted_Pi0] () {
            auto& node = fitted_Pi0->Get();
            // the returned bool is implicitly converted to double
            return ParticleTypeDatabase::Pi0.GetWindow(1).Contains(node.LVSum.M());
        });
    }


    // use mc_fake with complete 4pi (no lost photons)
    utils::MCFakeReconstructed mc_fake(true);

    unsigned nEvents = 0;
    unsigned nFailed = 0;

    while(true) {
        input::event_t event;
        if(!reader.ReadNextEvent(event))
            break;
        nEvents++;

        INFO("nEvents="+to_string(nEvents));

        TParticlePtr beam = event.MCTrue().ParticleTree->Get();
        auto mctrue_particles = mc_fake.Get(event.MCTrue());
        auto true_proton = mctrue_particles.Get(ParticleTypeDatabase::Proton).front();

        REQUIRE(mctrue_particles.GetAll().size() == 5);

        double prb = std_ext::NaN;
        TCandidatePtr best_proton;

        for(auto& p_proton : mctrue_particles.GetAll()) {
            auto& cand_proton = p_proton->Candidate;
            auto proton = make_shared<TParticle>(ParticleTypeDatabase::Proton, cand_proton);
            TParticleList photons;
            for(auto p_photon : mctrue_particles.GetAll()) {
                auto& cand_photon = p_photon->Candidate;
                if(cand_photon == cand_proton)
                    continue;
                photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, cand_photon));
            }

            // do the fit
            treefitter.PrepareFits(beam->Ek(), proton, photons);
            APLCON::Result_t res;
            unsigned nPerms = 0;
            while(treefitter.NextFit(res)) {
                nPerms++;
                if(res.Status != APLCON::Result_Status_t::Success)
                    continue;
                if(!std_ext::copy_if_greater(prb, res.Probability))
                    continue;
                best_proton = cand_proton;
            }
            if(true_proton->Candidate == cand_proton) {
                if(sort)
                    REQUIRE(nPerms==3);
                else
                    REQUIRE(nPerms==2);

            }
        }

        if(prb != Approx(1.0)) {
            nFailed++;
            continue;
        }

        REQUIRE(true_proton->Candidate == best_proton);
    }


    REQUIRE(nFailed == 3);
    REQUIRE(nEvents == 100);

}