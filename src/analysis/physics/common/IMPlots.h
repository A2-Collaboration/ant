#pragma once

#include "analysis/physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "utils/TriggerSimulation.h"
#include <vector>

class TH1D;
class TTree;

namespace ant {
namespace analysis {
namespace physics {

class IMPlots : public Physics {
public:
    utils::TriggerSimulation triggersimu;
    PromptRandom::Switch prs;
    std::vector<PromptRandom::Hist1> m;
    unsigned MinNGamma() const noexcept { return 2;}
    unsigned MaxNGamma() const noexcept { return unsigned(m.size())+2; }




public:
    IMPlots(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


class Symmetric2Gamma : public Physics {
protected:
    TH1D* h_symIM = nullptr;
    TTree* tree   = nullptr;

    double b_IM;
    double b_E;
    double b_E1;
    double b_E2;
    double b_theta1;
    double b_theta2;
    double b_phi1;
    double b_phi2;

    double perc = 0.2;

public:
    Symmetric2Gamma(const std::string& name, OptionsPtr opts);
    virtual ~Symmetric2Gamma();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

class IM_CB_TAPS_Plots : public Physics {
protected:

    struct hist_t {
        std::string prefix;

        TH1D* h_IM_All;
        TH1D* h_IM_CB;
        TH1D* h_IM_CB_ZeroVeto;
        TH1D* h_IM_CB_Pi0;
        TH1D* h_IM_TAPS;

        TH1D* h_Angle_CB;
        TH1D* h_Angle_TAPS;

        TH2D* h_ClusterHitTiming_CB;
        TH2D* h_ClusterHitTiming_TAPS;

        using range_t = interval<int>;
        const range_t n_CB;
        const range_t n_TAPS;
        static const range_t any;
        hist_t(const HistogramFactory& HistFac,
               const range_t& cb, const range_t& taps);
        void Fill(const TCandidatePtrList& c_CB, const TCandidatePtrList& c_TAPS) const;
        void ShowResult() const;


    };

    std::vector<hist_t> hists;

    TH1D* h_Mult_All;
    TH1D* h_Mult_CB;
    TH1D* h_Mult_TAPS;

public:
    IM_CB_TAPS_Plots(const std::string& name, OptionsPtr opts);
    virtual ~IM_CB_TAPS_Plots();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


}
}
}
