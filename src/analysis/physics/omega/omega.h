#pragma once

#include "analysis/physics/Physics.h"
#include "analysis/plot/SmartHist.h"
#include "analysis/utils/A2GeoAcceptance.h"

#include "base/interval.h"

#include <map>

class TH1D;
class TH2D;
class TH3D;

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

    static std::string GetDecayString(const data::ParticleList& particles);
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

    IntervalD omega_range = IntervalD(740,820);

    struct perDecayhists_t {
        TH1D* gg = nullptr;
        TH1D* ggg = nullptr;
        TH1D* mm = nullptr;
        TH1D* angle_p;
        TH1D* angle_p_ggg;
        TH1D* p_phi_diff;
        TH2D* calc_proton_energy_theta;
        TH2D* calc_proton_special;
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

}
}
}

std::string to_string(const ant::analysis::physics::OmegaBase::DataMode& m);
