#pragma once

#include "physics/Physics.h"
#include <map>
#include "base/WrapTTree.h"
#include "utils/TriggerSimulation.h"
#include "analysis/plot/PromptRandomHist.h"
#include "expconfig/detectors/Tagger.h"
#include "expconfig/detectors/TAPSVeto.h"
#include "expconfig/detectors/TAPS.h"
#include "calibration/converters/GeSiCa_SADC.h"
#include "calibration/converters/CATCH_TDC.h"

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

/**
 * @brief A playground class for checking taps stuff
 *
 */
class scratch_lheijken_checktaps: public Physics {
protected:

    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;

    std::shared_ptr<expconfig::detector::Tagger> tagger_detector;
    std::shared_ptr<expconfig::detector::TAPS> taps_detector;
    std::shared_ptr<expconfig::detector::TAPSVeto> veto_detector;

    TH1D* hTrigRefTiming;
    TH2D* hDRHUncalTimeAll;
    TH2D* hDRHCalTimeAll;
    TH2D* hDRHUncalTimeFirst;
    TH2D* hDRHCalTimeFirst;
    TH2D* hDRHUncalEnergy;
    TH2D* hDRHCalEnergy;
    TH2D* hDRHUncalTimeMult;
    TH2D* hDRHUncalEnergyMult;
    TH2D* hDRHUnCalEn_wTiming;
    TH2D* hDRHUnCalEn_woTiming;
    TH2D* hDRHCalEn_wTiming;
    TH2D* hDRHCalEn_woTiming;
    TH2D* hCHEnergy;
    TH2D* hCHTime;
    TH2D* hClTime;
    TH2D* hClEnergy;
    TH2D* hClTimeVsEnergy;
    TH2D* hCandClTime;
    TH2D* hCandClEnergy;
protected:
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);
    void CreateHistos();

public:
    scratch_lheijken_checktaps(const std::string& name, OptionsPtr opts);
    virtual ~scratch_lheijken_checktaps() {}

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;

};

}
}
}
