#pragma once

#include "physics/Physics.h"
#include "base/interval.h"
#include "base/piecewise_interval.h"
#include "plot/PromptRandomHist.h"
#include "utils/TriggerSimulation.h"

namespace ant {
namespace analysis {
namespace physics {

class MCEnergyThresholds : public Physics {

    const static BinSettings bins_Ek;
    const static BinSettings bins_clSize;
    const interval<double> CBThetaWindow;
    const PiecewiseInterval<double> CBHemisphereGap;

    utils::TriggerSimulation triggersimu;
    PromptRandom::Switch promptrandom;

    TH2D* h_nCaloClusters = nullptr;


public:

    class CBTAPS_t {
        const Detector_t::Type_t Type;
        const HistogramFactory HistFac;
        TH2D* h_Ecl_ClusterSize = nullptr;

    public:
        CBTAPS_t(Detector_t::Type_t type, const HistogramFactory& histFac);
        void Fill(const TCluster& caloCluster, double w) const;
        void Draw(canvas& c) const;
    };

private:

    const CBTAPS_t CB;
    const CBTAPS_t TAPS;

public:
    MCEnergyThresholds(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;

};


}}} // end of ant::analysis::physics
