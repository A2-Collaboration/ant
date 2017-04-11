#pragma once

#include "analysis/physics/Physics.h"
#include "base/piecewise_interval.h"
#include "base/WrapTTree.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {


struct ProcessTaggEff: public Physics {

    bool useTimeCut;

    unsigned seenEvents = 0;
    unsigned seenScalerBlocks = 0;
    unsigned nchannels = std::numeric_limits<unsigned>::quiet_NaN();

    TH1D* hist_scalers;
    TH1D* hist_tdchits;

    TH1D* hist_scalers_rate;
    TH1D* hist_tdchits_rate;

    TH1D* hist_tdc_times;
    TH2D* hist_tdc_times_ch;

    TH1D* hist_tdchits_wcut;
    TH1D* hist_tdc_times_wcut;
    TH2D* hist_tdc_times_ch_wcut;

    struct TreeScalarReads : WrapTTree {
        ADD_BRANCH_T(int,   nEvtsPerRead)
        ADD_BRANCH_T(TID,   EvID)

        ADD_BRANCH_T(double,    ExpLivetime)

        ADD_BRANCH_T(double,        Clock)
        ADD_BRANCH_T(double,        ExpTriggerRate)
        ADD_BRANCH_T(double,        PbRate)

        ADD_BRANCH_T(std::vector<int>,   TaggCounts)
        ADD_BRANCH_T(std::vector<double>,   TaggRates)

        ADD_BRANCH_T(std::vector<int>,   TDCCounts)
        ADD_BRANCH_T(std::vector<double>,   TDCRates)

        ADD_BRANCH_T(std::vector<std::vector<double>>,  TaggTimings)
    };

    TreeScalarReads scalerReads;

    static constexpr const char* treeName()        {return "scalerReads";}
    static constexpr const char* treeAccessName()  {return "ProcessTaggEff/scalerReads";}

    ProcessTaggEff(const std::string& name, OptionsPtr opts=nullptr);
    virtual ~ProcessTaggEff();

    virtual void ProcessEvent(const TEvent& ev, manager_t&) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

    void processBlock();
    void processTaggerHits(const TEvent& ev);
    void resetAll();
};

}}} // namespace ant::analysis::physics
