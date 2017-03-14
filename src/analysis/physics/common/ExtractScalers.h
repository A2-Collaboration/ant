#pragma once

#include "analysis/physics/Physics.h"
#include "base/WrapTTree.h"

#include "base/std_ext/math.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {


struct ExtractScalers: public Physics {

    using seconds_t = size_t;
    seconds_t   MediateOver;

    unsigned seenEvents = 0;
    unsigned seenScalerBlocks = 0;
    size_t nchannels;

    double clock=0;
    std::vector<double> TDCcounts;

    struct TreeScalers : WrapTTree {
        ADD_BRANCH_T(int,   nEvtsPerRead)

        ADD_BRANCH_T(double,    AbsTime)

        ADD_BRANCH_T(double,        IonChamber)
        ADD_BRANCH_T(double,        IonChamberError)
        ADD_BRANCH_T(double,        FaradayCup)
        ADD_BRANCH_T(double,        FaradayCupError)
        ADD_BRANCH_T(double,        EPTOr0)
        ADD_BRANCH_T(double,        EPTOr0Error)
        ADD_BRANCH_T(double,        EPTOr1)
        ADD_BRANCH_T(double,        EPTOr1Error)
        ADD_BRANCH_T(double,        EPTOr2)
        ADD_BRANCH_T(double,        EPTOr2Error)
        ADD_BRANCH_T(double,        EPTOr3)
        ADD_BRANCH_T(double,        EPTOr3Error)

        ADD_BRANCH_T(std::vector<double>,   TaggRates)
        ADD_BRANCH_T(std::vector<double>,   TaggRateErrors)
        ADD_BRANCH_T(std::vector<double>,   TDCRates)
        ADD_BRANCH_T(std::vector<double>,   TDCRateErrors)


    };
    struct medians_t {
        std_ext::RMS IonChamber;
        std_ext::RMS FaradayCup;
        std::vector<std_ext::RMS> EPTOrs;
        std::vector<std_ext::RMS> TaggRates;
        std::vector<std_ext::RMS> TDCRates;
        medians_t(const size_t nChannels):
            EPTOrs(4),
            TaggRates(nChannels),
            TDCRates(nChannels) {}
    };
    medians_t Medians;
    void calcRates();

    TreeScalers scalers;

    static constexpr const char* treeName()        {return "scalers";}
    static constexpr const char* treeAccessName()  {return "ExtractScalers/scalers";}
    static size_t nChannels();

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
