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

    SmartHistFactory HistFac_sig;
    SmartHistFactory HistFac_ref;

    TTree* treeSig;
    TTree* treeRef;

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

        sig_perDecayHists_t(SmartHistFactory& HistFac_parent, const std::string& decaystring);
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

        ref_perDecayHists_t(SmartHistFactory& HistFac_parent, const std::string& decaystring);
    };

    std::map<std::string, ref_perDecayHists_t> ref_perDecayHists;




    void ProcessSig(const data::ParticleTree_t& particletree, const data::Event::Data& data);
    void ProcessRef(const data::ParticleTree_t& particletree, const data::Event::Data& data);


public:
    EtapOmegaG(const std::string& name, PhysOptPtr opts);
    virtual void ProcessEvent(const data::Event& event) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};


}}}