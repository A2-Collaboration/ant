#ifndef OMEGA_H
#define OMEGA_H

#include "AntPhysics.h"
#include "base/interval.h"
#include "plot/SmartHist.h"
#include "A2GeoAcceptance.h"

class TH1D;
class TH2D;

namespace ant {
namespace analysis {

class Omega: public Physics {
protected:

    SmartHist1<const TLorentzVector&> omega_IM;
    SmartHist1<const TLorentzVector&> eta_IM;
    SmartHist1<const TLorentzVector&> p_MM;

    SmartHist1<int> omega_rec_multi;

    SmartHist1<int> nr_ngamma;
    SmartHist1<const TLorentzVector&> nr_3gim;
    SmartHist1<const TLorentzVector&> nr_2gim;

    int n=0;
    SmartHist1<const ParticlePtr&> test;

    IntervalD eta_im_cut;
    IntervalD pi0_im_cut;
    IntervalD omega_im_cut;
    IntervalD tagger_energy_cut;

    TLorentzVector target;

    SmartHist1<std::string> step_levels;

    SmartHist1< std::pair<const TLorentzVector&,const TLorentzVector&> > omega_mc_rec_angle;

    SmartHist1<const TLorentzVector&> makeInvMassPlot(const std::string& title, const std::string& xlabel, const std::string& ylabel, ant::BinSettings bins, const std::string& name="");
    SmartHist1< std::pair<const TLorentzVector&, const TLorentzVector&> > makeAngleDiffPlot(const std::string& title, const std::string& xlabel, const std::string& ylabel, BinSettings bins, const std::string& name);
    bool run_on_true;

public:
    Omega(const std::string& name, bool mctrue, const mev_t energy_scale=1000.0);
    virtual ~Omega() {}
    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};

class Omega2: public Physics {
protected:
    A2SimpleGeometry geo;

    double calcEnergySum(const ParticleList &particles) const;
    ParticleList getGeoAccepted(const ParticleList& p) const;
    SmartHist1<const TLorentzVector&> makeInvMassPlot(const std::string& title, const std::string& xlabel, const std::string& ylabel, ant::BinSettings bins, const std::string& name="");
    struct omega_decay {
        omega_decay(ParticlePtr omega_, ParticlePtr meson2_): omega(omega_), meson2(meson2_), score(0.0) {}
        ParticlePtr omega;
        ParticlePtr meson2;
        double score;
    };

    struct settings_t {
        double esum_threshold;
        IntervalD omega_IM_cut;
        IntervalD eta_IM_cut;
        IntervalD pi0_IM_cut;
    };

    settings_t settings;

    bool run_on_true;

    SmartHist1<std::string> step_levels;

    SmartHist1<const TLorentzVector&> omega_IM;
    SmartHist1<const TLorentzVector&> eta_IM;
    SmartHist1<const TLorentzVector&> p_MM;
    SmartHist1<const TLorentzVector&> pi0_IM;
    SmartHist1<const TLorentzVector&> gggIM;
    SmartHist1<const TLorentzVector&>  ggIM;

    SmartHist1<std::string>  found_candidates;

public:
    Omega2(const std::string& name, bool mctrue);
    virtual ~Omega2() {}
    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};

}
}
#endif
