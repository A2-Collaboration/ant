#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/plot/SmartHist.h"
#include "analysis/utils/A2GeoAcceptance.h"
#include "base/Tree.h"
#include "base/interval.h"
#include "analysis/utils/particle_tools.h"
#include "base/std_ext/math.h"

#include "analysis/utils/KinFitter.h"
#include "base/interval.h"
#include "analysis/plot/PromptRandomHist.h"

#include <map>

class TH1D;
class TH2D;
class TH3D;
class TTree;
namespace ant {

namespace analysis {
namespace physics {

class OmegaMCTruePlots: public Physics {
public:
    struct PerChannel_t {
        std::string title;
        TH2D* proton_E_theta = nullptr;

        PerChannel_t(const std::string& Title, SmartHistFactory& hf);

        void Show();
        void Fill(const data::Event::Data& d);
    };

    std::map<std::string,PerChannel_t> channels;

    OmegaMCTruePlots(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event& event);
    void Finish();
    void ShowResult();
};

class OmegaBase: public Physics {

public:
    enum class DataMode {
        MCTrue,
        Reconstructed
    };

protected:
    utils::A2SimpleGeometry geo;
    double calcEnergySum(const data::ParticleList &particles) const;
    data::ParticleList getGeoAccepted(const data::ParticleList& p) const;

    DataMode mode = DataMode::Reconstructed;

    virtual void Analyse(const data::Event::Data& data, const data::Event& event) =0;



public:
    OmegaBase(const std::string &name, PhysOptPtr opts);
    virtual ~OmegaBase() = default;

    virtual void ProcessEvent(const data::Event& event) override;
    void Finish() override;
    void ShowResult() override;


};

class OmegaEtaG: public OmegaBase {

protected:

    TH2D* ggg_gg;
    TH2D* ggg_gg_bg;    // if not from omega decay
    TH2D* ggg_gg_all;
    //TH1D* gg_bg;        // if from omega but not pi0 or eta decay

    TH1D* ggg;
    TH1D* ggg_omega;
    TH1D* ggg_bg;
    TH1D* ggg_omega_pi0oreta;

    TH2D* ggg_gg_omega_eta;
    TH2D* ggg_gg_omega_pi0;

    TH1D* steps;

    IntervalD omega_range = IntervalD(680,780);

    struct perDecayhists_t {
        TH1D* gg = nullptr;
        TH1D* ggg = nullptr;
        TH1D* mm = nullptr;
        TH1D* angle_p;
        TH1D* angle_p_ggg;
        TH1D* p_phi_diff;
        TH2D* calc_proton_energy_theta;
        TH2D* calc_proton_special;
        TH1D* nCand;
    };

    perDecayhists_t makePerDecayHists(const std::string &title="");

    std::map<std::string, perDecayhists_t> gg_decays;

    virtual void Analyse(const data::Event::Data& data, const data::Event& event) override;

    BinSettings imbinning = BinSettings(1000);
    BinSettings mmbinning = BinSettings(1000, 400,1400);

public:
    OmegaEtaG(const std::string& name, PhysOptPtr opts);
    virtual ~OmegaEtaG() = default;
    void ShowResult() override;
};



class OmegaMCTree : public Physics {
protected:
    TTree* tree = nullptr;
    TLorentzVector proton_vector;
    TLorentzVector omega_vector;
    TLorentzVector eta_vector;
    TLorentzVector gamma1_vector;
    TLorentzVector gamma2_vector;
    TLorentzVector gamma3_vector;

public:

    OmegaMCTree(const std::string& name, PhysOptPtr opts);
    virtual ~OmegaMCTree();

    void ProcessEvent(const data::Event& event) override;
    void ShowResult() override;
    TLorentzVector getGamma1() const;
    void setGamma1(const TLorentzVector& value);
};


class OmegaEtaG2 : public OmegaBase {


    // OmegaBase interface
protected:
    void Analyse(const data::Event::Data &data, const data::Event &event) override;


    enum SigBgFlag_t {
        flagSignal = 0,
        flagReference,
        flagBackground
    };

    SigBgFlag_t identify(const data::Event &event) const;


    std::shared_ptr<ant::Tree<const ParticleTypeDatabase::Type&>> signal_tree;
    std::shared_ptr<ant::Tree<const ParticleTypeDatabase::Type&>> reference_tree;


    //======== Tree ===============================================================

    TTree*  tree = nullptr;

    analysis::utils::ParticleVars b_g1;
    analysis::utils::ParticleVars b_g2;
    analysis::utils::ParticleVars b_g3;
    analysis::utils::ParticleVars b_p;
    analysis::utils::ParticleVars b_ggg;
    analysis::utils::ParticleVars b_mmvector;

    double b_pTime   = {};
    double b_gggTime = {};
    double b_ggIM[3] = {};

    double b_copl_angle = 0.0;
    double b_p_mm_angle = 0.0;

    int    b_TagCh     = -1;
    double b_TagTime   = 0.0;
    double b_TagE      = 0.0;
    double b_TagW      = 0.0;


    int       b_SigBgFlag = flagSignal;

    double b_ggIM_real    = {};
    double b_ggIM_comb[2] = {};

    double b_BachelorE[3] = {};

    double b_CBAvgTime = 0.0;


    double kinfit_chi2       = 0.0;
    bool   b_fitok           = false;
    unsigned b_fitIterations = 0;


    //======== Settings ===========================================================

    bool data_proton = true;
    bool data_tagger = true;

    double cut_ESum = 550.0;
    double cut_Copl = std_ext::degree_to_radian(15.0);

    interval<double> photon_E_cb   = { 50.0,1600.0};
    interval<double> photon_E_taps = {200.0,1600.0};
    interval<double> proton_theta  = std_ext::degree_to_radian(interval<double>({0.0, 45.0}));

    double calcEnergySum2(const data::Event::Data &e) const;

    struct expected_peak_t {
        double Mean;
        double Sigma;
        expected_peak_t(double mean, double sigma) :
            Mean(mean), Sigma(sigma) {}
    };

    const expected_peak_t omega_peak = {7.64266e+02, 3.29420e+01};
    const expected_peak_t eta_peak   = {5.30815e+02, 2.93928e+01};
    const expected_peak_t pi0_peak   = {1.31331e+02, 1.04835e+01};


    TH1D* steps;

    TH1D* h_TotalEvents;

    std::map<std::string, TH1D*> pulls;

    const std::map<short, std::string> component = {{0, "Energy"}, {1, "Theta"}, {2, "Phi"}};
    const BinSettings pull_bins = BinSettings(201, -10, 10);


    ant::analysis::PromptRandom::Switch promptrandom;
    utils::KinFitter fitter;

    bool AcceptedPhoton(const data::ParticlePtr& photon);
    bool AcceptedProton(const data::ParticlePtr& proton);

    data::ParticleList FilterPhotons(const data::ParticleList& list);
    data::ParticleList FilterProtons(const data::ParticleList& list);

public:
    OmegaEtaG2(const std::string& name, PhysOptPtr opts);
    virtual ~OmegaEtaG2();

};

}
}
}

std::string to_string(const ant::analysis::physics::OmegaBase::DataMode& m);
