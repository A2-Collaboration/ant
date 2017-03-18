#pragma once

#include "physics/Physics.h"
#include "plot/PromptRandomHist.h"
#include "utils/fitter/KinFitter.h"
#include "utils/Uncertainties.h"
#include "utils/TriggerSimulation.h"

namespace ant {
namespace analysis {
namespace physics {

class TriggerSimulation : public Physics {
protected:
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;

    TH1D* steps;

    TH1D* h_CBESum_raw;
    TH1D* h_CBESum_pr;
    TH1D* h_CBTiming;
    TH2D* h_CBTiming_CaloE;

    struct ClusterPlots_t {
        ClusterPlots_t(const HistogramFactory& HistFac);
        void Fill(const TEventData& recon) const;
        void Show(canvas& c) const;
        TH2D* h_CaloE_ClSize;
        TH2D* h_CaloE_nCl;
        TH2D* h_CaloE_Time;
        TH1D* h_Hits_stat;
        TH2D* h_Hits_E_t;
    };

    const ClusterPlots_t Clusters_All;
    const ClusterPlots_t Clusters_Tail;


    TH1D* h_TaggT;
    TH1D* h_zvertex;


    utils::UncertaintyModelPtr fit_model;
    utils::KinFitter fitter;


public:
    TriggerSimulation(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};


}}}