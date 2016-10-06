#pragma once

#include "analysis/physics/Physics.h"
#include "base/WrapTTree.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {


struct ExtractScalers: public Physics {

    using seconds_t = size_t;
    seconds_t   MediateOver;

    unsigned seenEvents = 0;
    unsigned seenScalerBlocks = 0;
    unsigned nchannels = std::numeric_limits<unsigned>::quiet_NaN();

    const HistogramFactory histFac;

    TH1D* hist_scalers;
    TH1D* hist_tdchits;

    TH1D* hist_scalers_rate;
    TH1D* hist_tdchits_rate;

    double clock=0;
    std::vector<double> TDCcounts;

    struct TreeScalers : WrapTTree {
        ADD_BRANCH_T(int,   nEvtsPerRead)

        ADD_BRANCH_T(double,    AbsTime)

        ADD_BRANCH_T(double,        IonChamber)
        ADD_BRANCH_T(double,        FaradayCup)
        ADD_BRANCH_T(double,        EPTOr0)
        ADD_BRANCH_T(double,        EPTOr1)
        ADD_BRANCH_T(double,        EPTOr2)
        ADD_BRANCH_T(double,        EPTOr3)

        ADD_BRANCH_T(std::vector<double>,   TaggRates)
        ADD_BRANCH_T(std::vector<double>,   TDCRates)


    };

    TreeScalers scalers;

    static constexpr auto treeName()        {return "scalers";}
    static constexpr auto treeAccessName()  {return "ExtractScalers/scalers";}

    ExtractScalers(const std::string& name, OptionsPtr opts=nullptr);
    virtual ~ExtractScalers();

    virtual void ProcessEvent(const TEvent& ev, manager_t&) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

    void processBlock();
    void processTaggerHits(const TEvent& ev);
    void resetAll();
};

}}} // namespace ant::analysis::physics
