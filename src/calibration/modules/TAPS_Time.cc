#include "TAPS_Time.h"

#include "calibration/fitfunctions/FitGausPol0.h"

#include "expconfig/detectors/TAPS.h"

#include "base/Logger.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;

TAPS_Time::TAPS_Time(shared_ptr<expconfig::detector::TAPS> taps,
           shared_ptr<DataManager> calmgr,
           Calibration::Converter::ptr_t converter,
           const interval<double>& timeWindow_BaF2, // default {-inf, inf}
           const interval<double>& timeWindow_PbWO4 // default {-inf, inf}
           ) :
    Time(taps,
         calmgr,
         converter,
         -170, // for BaF2
         std::make_shared<calibration::gui::FitGausPol0>(),
         timeWindow_BaF2, // for BaF2
         -0.100 // for BaF2
         )
{
    // DefaultGains, DefaultOffsets and TimeWindows are different for
    // PbWO4 elements

    for(unsigned ch=0; ch<taps->GetNChannels(); ch++)
    {
        if(taps->IsPbWO4(ch)) {
            DefaultOffsets[ch] = 400;
            DefaultGains[ch] = 0.1;
            TimeWindows[ch] = timeWindow_PbWO4;
        }
    }
}

std::unique_ptr<Physics> TAPS_Time::GetPhysicsModule()
{
    return std_ext::make_unique<ThePhysics>(GetName(), "Offsets", Detector);
}


TAPS_Time::ThePhysics::ThePhysics(const string& name,
                                  const string& histName,
                                  const std::shared_ptr<Detector_t>& theDetector) :
    Time::ThePhysics(name, histName, theDetector)
{
    string detectorName(Detector_t::ToString(detector->Type));
    hTimeToTagger = HistFac.makeTH2D(
                        detectorName + string(" - Time relative to tagger"),
                        "time [ns]",
                        detectorName + " channel",
                        BinSettings(2000,-1500,1500),
                        BinSettings(detector->GetNChannels()),
                        "hTimeToTagger"
                        );
}

void TAPS_Time::ThePhysics::ProcessEvent(const Event& event)
{
    Time::ThePhysics::ProcessEvent(event);
    for(const auto& cand : event.Reconstructed().Candidates()) {
        for(const auto& cluster: cand->Clusters) {
            if(cluster.Detector != detector->Type)
                continue;
            for(const auto& taggerhit : event.Reconstructed().TaggerHits()) {
                const double relative_time = cluster.Time - taggerhit->Time();
                hTimeToTagger->Fill(relative_time, cluster.CentralElement);
            }
        }
    }
}

void TAPS_Time::ThePhysics::ShowResult()
{
    Time::ThePhysics::ShowResult();
    canvas(GetName()+" to Tagger")  << drawoption("colz") << hTimeToTagger << endc;
}
