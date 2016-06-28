#include "PID_Energy.h"

#include "gui/CalCanvas.h"

#include "calibration/fitfunctions/FitGaus.h"

#include "expconfig/detectors/PID.h"

#include "base/Logger.h"

#include "tree/TEventData.h"

#include <list>


using namespace std;
using namespace ant;
using namespace ant::calibration;

template <typename Detector>
bool checkDefaultArraySize(Detector det, const std::vector<double> arr) {
    return (arr.size() == 1 || arr.size() == det->GetNChannels());
}

PID_Energy::PID_Energy(
        std::shared_ptr<expconfig::detector::PID> pid,
        std::shared_ptr<DataManager> calmgr,
        Calibration::Converter::ptr_t converter,
        const std::vector<double>& defaultPedestals,
        const std::vector<double>& defaultGains,
        const std::vector<double>& defaultThresholds,
        const std::vector<double>& defaultRelativeGains
        ) :
    Energy(pid->Type,
           calmgr,
           converter,
           defaultPedestals,
           defaultGains,
           defaultThresholds,
           defaultRelativeGains),
    pid_detector(pid)
{
    if(!checkDefaultArraySize(pid, defaultPedestals)) {
        throw std::runtime_error("PID default pedestals size wrong");
    }
    if(!checkDefaultArraySize(pid, defaultGains)) {
        throw std::runtime_error("PID default gains size wrong");
    }
    if(!checkDefaultArraySize(pid, defaultThresholds)) {
        throw std::runtime_error("PID default threshols size wrong");
    }
    if(!checkDefaultArraySize(pid, defaultRelativeGains)) {
        throw std::runtime_error("PID default relative size wrong");
    }
}


PID_Energy::~PID_Energy()
{

}

void PID_Energy::GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >& guis)
{
    guis.emplace_back(std_ext::make_unique<GUI_Pedestals>(
                          GetName(),
                          Pedestals,
                          calibrationManager,
                          pid_detector,
                          make_shared<gui::FitGaus>()
                          ));
    guis.emplace_back(std_ext::make_unique<GUI_Banana>(
                          GetName(),
                          RelativeGains,
                          calibrationManager,
                          pid_detector
                          ));

}

void ant::calibration::PID_Energy::ApplyTo(TEventData& reconstructed)
{
    // search for PID/CB candidates and correct PID energy by CB theta angle
    for(TCandidate& cand : reconstructed.Candidates) {
        const bool cb_and_pid = cand.Detector & Detector_t::Type_t::CB &&
                                cand.Detector & Detector_t::Type_t::PID;
        if(!cb_and_pid)
            continue;
        cand.VetoEnergy *= sin(cand.Theta);
    }
}
