#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/plot/SmartHist.h"
#include "analysis/utils/A2GeoAcceptance.h"
#include "base/Tree.h"
#include "base/interval.h"
#include "analysis/utils/particle_tools.h"
#include "base/std_ext/math.h"

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

    struct mParticleVars : ant::analysis::utils::ParticleVars {
        double matchAngle = -1.0;

        virtual void SetBranches(TTree* tree, const std::string& name) override;

        using ParticleVars::ParticleVars;

        virtual ~mParticleVars() {}
    };

    std::shared_ptr<ant::Tree<const ParticleTypeDatabase::Type&>> signal_tree;
    std::shared_ptr<ant::Tree<const ParticleTypeDatabase::Type&>> reference_tree;

    analysis::utils::ParticleVars pbranch;
    double pTime = {};

    analysis::utils::ParticleVars gggbranch;
    double gggTime = {};

    double ggIM[3] = {};

    analysis::utils::ParticleVars calcp;

    double angle_p_calcp = {};

    mParticleVars g1branch;
    mParticleVars g2branch;
    mParticleVars g3branch;

    int    tagch = -1;
    double tagtime = {};

    int    rf = -1;

    TTree*  tree = nullptr;

    bool data_proton = true;
    bool data_tagger = true;
    double ESum_cut = 550.0;

    double calcEnergySum2(const data::Event::Data &e) const;

    enum class channel_type_t {
        Signal,
        Reference,
        Background
    };

    channel_type_t identify(const data::Event& event) const;



    struct particleCuts_t {
        interval<double> E_range = {0,1600};
        interval<double> Theta_range = {0, 2*M_PI};

        bool TestParticle(const data::Particle& p) const;
    };

    particleCuts_t photon_cut;
    particleCuts_t proton_cut;

    data::ParticleList FilterParticles(const data::ParticleList& list, const particleCuts_t& cuts) const;

    struct expected_peak_t {
        double Mean;
        double Sigma;
        expected_peak_t(double mean, double sigma) :
            Mean(mean), Sigma(sigma) {}
    };

    const expected_peak_t omega_peak = {7.64266e+02, 3.29420e+01};
    const expected_peak_t eta_peak   = {5.30815e+02, 2.93928e+01};
    const expected_peak_t pi0_peak   = {1.31331e+02, 1.04835e+01};

    double Chi2_Omega = 0.0;
    double Chi2_Pi0[3] = {};
    double Chi2_Eta[3] = {};
    char   bestEtaIn = 0;
    char   bestPi0In = 0;

    double   bestChi = 0;
    char     bestHyp   = 0; // 1== Eta, 2==Pi0, 0==??

    double EgOmegaSys[3] = {};

    TH1D* steps;

    const interval<double> complcut = std_ext::degree_to_radian(interval<double>::CenterWidth(180,30));

    TH1D* h_TotalEvents;


public:
    OmegaEtaG2(const std::string& name, PhysOptPtr opts);
    virtual ~OmegaEtaG2();

};

}
}
}

std::string to_string(const ant::analysis::physics::OmegaBase::DataMode& m);
