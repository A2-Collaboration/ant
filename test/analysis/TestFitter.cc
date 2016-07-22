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
        utils::Uncertainties_t  u{
                    0.05*particle.Ek(),
                    std_ext::degree_to_radian(2.0),
                    std_ext::degree_to_radian(2.0)
        };
        if(ProtonUnmeasured && particle.Type() == ParticleTypeDatabase::Proton)
            u.sigmaE = 0;
        return u;
    }
};

struct RMS_t {
    unsigned n = 0;
    double sum = 0;
    double sum2 = 0;
    void Add(double v) {
        ++n;
        sum += v;
        sum2 += std_ext::sqr(v);
    }
    double GetMean() const {
        return sum/n;
    }
    double GetRMS() const {
        return std::sqrt( sum2/n - std_ext::sqr(GetMean()) );
    }
};

struct Pulls_t {
    RMS_t Ek;
    RMS_t Theta;
    RMS_t Phi;

    void Fill(const utils::Fitter::FitParticle& p) {
        Ek.Add(p.Ek.Pull);
        Theta.Add(p.Theta.Pull);
        Phi.Add(p.Phi.Pull);
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
        kinfitter.SetZVertexSigma(0.0);
        test::EnsureSetup(); // needed for MCFake
    }

    // use mc_fake with complete 4pi (no lost photons)
    auto mc_fake = z_vertex ? std_ext::make_unique<utils::MCFakeReconstructed>(true) : nullptr;
    auto mc_smear = smeared ? std_ext::make_unique<utils::MCSmear>(model) : nullptr;

    unsigned nEvents = 0;
    unsigned nFitOk = 0;
    unsigned nFitIterations = 0;

    RMS_t pulls_Beam;
    Pulls_t pulls_Photons;
    Pulls_t pulls_Proton;

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
            REQUIRE(kinfitter.GetFittedZVertex() == Approx(0.0).epsilon(6e-3));
        else
            REQUIRE(std::isnan(kinfitter.GetFittedZVertex()));

        const auto& fitparticles = kinfitter.GetFitParticles();
        REQUIRE(fitparticles.size() == 3);
        REQUIRE(fitparticles.front().Particle->Type() == ParticleTypeDatabase::Proton);
        auto it_fitparticle = fitparticles.begin();
        pulls_Proton.Fill(*it_fitparticle);
        ++it_fitparticle;
        while (it_fitparticle != fitparticles.end()) {
            pulls_Photons.Fill(*it_fitparticle);
            ++it_fitparticle;
        }
    }

    CHECK(nEvents==1000);
    CHECK(nFitOk==nEvents);

    if(smeared) {
        CHECK(nFitIterations > nEvents*2); // fitter should take longer to converge
        if(proton_unmeas) {
            CHECK(pulls_Proton.Ek.GetMean() == Approx(0));
            CHECK(pulls_Proton.Ek.GetRMS() == Approx(0));
        }
        else {
            CHECK(pulls_Proton.Ek.GetMean() == Approx(0).epsilon(0.06));
            CHECK(pulls_Proton.Ek.GetRMS() == Approx(1).epsilon(0.15));
        }

        CHECK(pulls_Proton.Theta.GetMean() == Approx(0).epsilon(0.06));
        CHECK(pulls_Proton.Theta.GetRMS() == Approx(1).epsilon(0.06));
        CHECK(pulls_Proton.Phi.GetMean() == Approx(0).epsilon(0.04));
        CHECK(pulls_Proton.Phi.GetRMS() == Approx(1).epsilon(0.02));

        CHECK(pulls_Photons.Ek.GetMean() == Approx(0).epsilon(0.04));
        CHECK(pulls_Photons.Ek.GetRMS() == Approx(1).epsilon(0.01));
        CHECK(pulls_Photons.Theta.GetMean() == Approx(0).epsilon(0.03));
        CHECK(pulls_Photons.Theta.GetRMS() == Approx(1).epsilon(0.04));
        CHECK(pulls_Photons.Phi.GetMean() == Approx(0).epsilon(0.003));
        CHECK(pulls_Photons.Phi.GetRMS() == Approx(1).epsilon(0.01));
    }
    else {
        // unsmeared, so all pulls should be delta peaks...
        constexpr double eps = 1e-3;
        CHECK(pulls_Proton.Ek.GetMean() == Approx(0).epsilon(eps));
        CHECK(pulls_Proton.Ek.GetRMS() == Approx(0).epsilon(eps));
        CHECK(pulls_Proton.Theta.GetMean() == Approx(0).epsilon(eps));
        CHECK(pulls_Proton.Theta.GetRMS() == Approx(0).epsilon(eps));
        CHECK(pulls_Proton.Phi.GetMean() == Approx(0).epsilon(eps));
        CHECK(pulls_Proton.Phi.GetRMS() == Approx(0).epsilon(eps));

        CHECK(pulls_Photons.Ek.GetMean() == Approx(0).epsilon(eps));
        CHECK(pulls_Photons.Ek.GetRMS() == Approx(0).epsilon(eps));
        CHECK(pulls_Photons.Theta.GetMean() == Approx(0).epsilon(eps));
        CHECK(pulls_Photons.Theta.GetRMS() == Approx(0).epsilon(eps));
        CHECK(pulls_Photons.Phi.GetMean() == Approx(0).epsilon(eps));
        CHECK(pulls_Photons.Phi.GetRMS() == Approx(0).epsilon(eps));
    }
}
