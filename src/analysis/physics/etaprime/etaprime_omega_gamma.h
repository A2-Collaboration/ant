#pragma once

#include "analysis/physics/Physics.h"
#include "utils/particle_tools.h"

class TH1D;
class TH2D;
class TH3D;

namespace ant {
namespace analysis {
namespace physics {

class EtapOmegaG : public Physics {

    struct expected_peak_t {
        double Mean;
        double Sigma;
        expected_peak_t(double mean, double sigma) :
            Mean(mean), Sigma(sigma) {}
    };

    // means/sigma extracted from gg/ggg/gggg histograms for signal channel
    const expected_peak_t Pi0 = {126, 15};
    const expected_peak_t Omega = {735, 32};
    const expected_peak_t EtaPrime_sig = {895, 27};
    // extracted from gg histogram for reference channel
    const expected_peak_t EtaPrime_ref = {905, 29};

    TH1D* sig_steps;
    TH1D* ref_steps;


    struct sig_perDecayHists_t {
        TH1D* gggg;
        TH1D* ggg;
        TH1D* gg;

        TH1D* Proton_Copl;

        TH2D* IM_etap_omega;
        TH1D* IM_pi0;

        TH1D* MM_gggg;
        TH1D* MM_etap;

        TH2D* Proton_ThetaPhi;
        TH1D* Proton_Energy;

        TH1D* Chi2_All;
        TH1D* Chi2_Best;

        sig_perDecayHists_t(SmartHistFactory& HistFac_parent, const std::string& decaystring,
                            const std::string& prefix);
    };

    std::map<std::string, sig_perDecayHists_t> sig_perDecayHists;

    struct ref_perDecayHists_t {
        TH1D* gg;

        TH2D* Proton_ThetaPhi;
        TH1D* Proton_Energy;

        TH1D* IM_etap;

        TH1D* MM_etap;

        TH1D* Proton_Copl;

        bool IsRelevant() const {
            return IM_etap->GetEntries()>0;
        }

        ref_perDecayHists_t(SmartHistFactory& HistFac_parent, const std::string& decaystring,
                            const std::string& prefix);
    };

    std::map<std::string, ref_perDecayHists_t> ref_perDecayHists;


    template<typename T>
    T& getHistogram(const std::string& prefix,
                    const data::ParticleTree_t& particletree, std::map<std::string, T>& perDecayHists) {
        const std::string& decaystring = utils::ParticleTools::GetDecayString(particletree);
        // search map only once even on insert
        auto it_h = perDecayHists.lower_bound(decaystring);
        if(it_h == perDecayHists.end() || perDecayHists.key_comp()(decaystring, it_h->first)) {
            // decaystring does not exist
            it_h = perDecayHists.emplace_hint(it_h, decaystring, T(HistFac, decaystring, prefix));
        }
        return it_h->second;
    }

    void ProcessSig(const data::ParticleTree_t& particletree, const data::Event::Data& data);
    void ProcessRef(const data::ParticleTree_t& particletree, const data::Event::Data& data);


public:
    EtapOmegaG(PhysOptPtr opts);
    virtual void ProcessEvent(const data::Event& event) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};


}}}