#include "Time.h"

#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::data;
using namespace ant::analysis::physics;

Time::Time(const Detector_t::Type_t& detectorType,
           const string& name, analysis::PhysOptPtr opts)
    :
      Physics(name, opts),
      DetectorType(detectorType)
{
    auto detector = ExpConfig::Setup::GetDetector(DetectorType);
    string detectorName(Detector_t::ToString(DetectorType));
    hTime = HistFac.makeTH2D(detectorName + string(" - Time"),
                             "time [ns]",
                             detectorName + " channel",
                             BinSettings(2000,-400,400),
                             BinSettings(detector->GetNChannels()),
                             "Time"
                             );
    hTimeToTagger = HistFac.makeTH2D(
                        detectorName + string(" - Time relative to tagger"),
                        "time [ns]",
                        detectorName + " channel",
                        BinSettings(2000,-1500,1500),
                        BinSettings(detector->GetNChannels()),
                        "hTimeToTagger"
                        );

    // handle tagger differently
    isTagger = dynamic_pointer_cast<TaggerDetector_t, Detector_t>(detector) != nullptr;
}

void Time::ProcessEvent(const Event& event)
{
    // handle Tagger differently
    if(isTagger)
    {
        for (const auto& tHit: event.Reconstructed.TaggerHits)
            hTime->Fill(tHit.Time,tHit.Channel);
        return;
    }

    for(const auto& cand: event.Reconstructed.Candidates) {
        for(const auto& cluster: cand->Clusters) {
            if(cluster.Detector != DetectorType)
                continue;
            hTime->Fill(cluster.Time,cluster.CentralElement);
            for(const auto& taggerhit : event.Reconstructed.TaggerHits) {
                const double relative_time = cluster.Time - taggerhit.Time;
                hTimeToTagger->Fill(relative_time, cluster.CentralElement);
            }
        }
    }
}

void Time::ShowResult()
{
    canvas(GetName())
            << drawoption("colz") << hTime
            << drawoption("colz") << hTimeToTagger
            << endc;
}

namespace ant {
namespace analysis {
namespace physics {

struct EPT_Time : Time {
    EPT_Time(const std::string& name, PhysOptPtr opts) :
        Time(Detector_t::Type_t::EPT, name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(EPT_Time)

struct CB_Time : Time {
    CB_Time(const std::string& name, PhysOptPtr opts) :
        Time(Detector_t::Type_t::CB, name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(CB_Time)

struct PID_Time : Time {
    PID_Time(const std::string& name, PhysOptPtr opts) :
        Time(Detector_t::Type_t::PID, name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(PID_Time)

struct TAPS_Time : Time {
    TAPS_Time(const std::string& name, PhysOptPtr opts) :
        Time(Detector_t::Type_t::TAPS, name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(TAPS_Time)

struct TAPSVeto_Time : Time {
    TAPSVeto_Time(const std::string& name, PhysOptPtr opts) :
        Time(Detector_t::Type_t::TAPSVeto, name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(TAPSVeto_Time)


}}}








