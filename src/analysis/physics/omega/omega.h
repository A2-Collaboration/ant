#ifndef OMEGA_H
#define OMEGA_H

#include "physics/Physics.h"
#include "base/interval.h"
#include "plot/SmartHist.h"
#include "A2GeoAcceptance.h"
#include <map>

class TH1D;
class TH2D;
class TH3D;

namespace ant {
namespace analysis {

class OmegaBase: public Physics {

public:
    enum class DataMode {
        MCTrue,
        Reconstructed
    };

protected:
    A2SimpleGeometry geo;
    double calcEnergySum(const ParticleList &particles) const;
    ParticleList getGeoAccepted(const ParticleList& p) const;

    DataMode mode = DataMode::Reconstructed;

    virtual void Analyse(const Event::Data& data, const Event& event) =0;

    static std::string GetDecayString(const ParticleList& particles);

public:
    OmegaBase(const string &name, const DataMode m);
    virtual ~OmegaBase() = default;

    virtual void ProcessEvent(const Event& event) override;
    void Finish() override;
    void ShowResult() override;

};

class OmegaEtaG: public OmegaBase {

protected:

    TH2D* ggg_gg;
    TH3D* ggg_gg_decays;
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
    std::map<std::string, TH1D*> gg_decays;

    virtual void Analyse(const Event::Data& data, const Event& event) override;

    BinSettings imbinning = BinSettings(1000);

public:
    OmegaEtaG(DataMode m);
    virtual ~OmegaEtaG() = default;
    void ShowResult() override;
};

}
}

std::string to_string(const ant::analysis::OmegaBase::DataMode& m);
#endif
