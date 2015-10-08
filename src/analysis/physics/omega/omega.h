#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/plot/SmartHist.h"
#include "analysis/utils/A2GeoAcceptance.h"

#include "base/interval.h"

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

    OmegaMCTruePlots(PhysOptPtr opts);

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
    OmegaEtaG(PhysOptPtr opts);
    virtual ~OmegaEtaG() = default;
    void ShowResult() override;
};



class OmegeMCTree : public Physics {
protected:
    TTree* tree = nullptr;
    TLorentzVector proton_vector;
    TLorentzVector omega_vector;
    TLorentzVector eta_vector;
    TLorentzVector gamma1_vector;
    TLorentzVector gamma2_vector;
    TLorentzVector gamma3_vector;

public:

    OmegeMCTree(PhysOptPtr opts);
    virtual ~OmegeMCTree();

    void ProcessEvent(const data::Event& event) override;
    void ShowResult() override;
    TLorentzVector getGamma1() const;
    void setGamma1(const TLorentzVector& value);
};


class OmegaEtaG2 : public OmegaBase {


    // OmegaBase interface
protected:
    void Analyse(const data::Event::Data &data, const data::Event &event) override;

    double pEk = {};
    double pTheta = {};
    double pPhi = {};
    double pTime = {};

    double gggIM = {};
    double gggTheta = {};
    double gggPhi = {};
    double gggTime = {};
    double gggE = {};

    double ggIM[3] = {};
    double MM = {};

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

public:
    OmegaEtaG2(PhysOptPtr opts);
    virtual ~OmegaEtaG2();

};

}
}
}

std::string to_string(const ant::analysis::physics::OmegaBase::DataMode& m);
