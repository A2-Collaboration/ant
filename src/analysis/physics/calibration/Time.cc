#include "Time.h"

#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

Time::Time(const Detector_t::Type_t& detectorType,
           const string& name, OptionsPtr opts)
    :
      Physics(name, opts)
{
    Detector = ExpConfig::Setup::GetDetector(detectorType);
    string detectorName(Detector_t::ToString(Detector->Type));
    hTime = HistFac.makeTH2D(detectorName + " - Time",
                             "time [ns]",
                             detectorName + " channel",
                             BinSettings(2000,-400,400),
                             BinSettings(Detector->GetNChannels()),
                             "Time"
                             );
    hTimeToF = HistFac.makeTH2D(detectorName + " - Time for ToF",
                             "time [ns]",
                             detectorName + " channel",
                             BinSettings(1000,-50,50),
                             BinSettings(Detector->GetNChannels()),
                             "Time_ToF" // for TAPS_ToF offsets...
                             );
    hTimeToTagger = HistFac.makeTH2D(
                        detectorName + " - Time relative to tagger",
                        "time [ns]",
                        detectorName + " channel",
                        BinSettings(2000,-1500,1500),
                        BinSettings(Detector->GetNChannels()),
                        "hTimeToTagger"
                        );
    hCBTriggerTiming = HistFac.makeTH1D("CB - Energy-averaged time",
                                "time [ns]",
                                "#",
                                BinSettings(500,-15,15),
                                "hCBTriggerTiming");
    hTimeMultiplicity = HistFac.makeTH2D(detectorName + " - Time Hit Multiplicity",
                                         "multiplicity",
                                         detectorName + " channel",
                                         BinSettings(8),
                                         BinSettings(Detector->GetNChannels()),
                                         "hTimeMultiplicity"
                                         );


    // handle tagger differently
    isTagger = dynamic_pointer_cast<TaggerDetector_t, Detector_t>(Detector) != nullptr;

}

void Time::ProcessEvent(const TEvent& event, manager_t&)
{
    const double CBTimeAvg = event.Reconstructed().Trigger.CBTiming;
    hCBTriggerTiming->Fill(CBTimeAvg);

    std::map<unsigned, unsigned> multiplicity;

    // handle Tagger differently
    if(isTagger)
    {
        for (const auto& tHit: event.Reconstructed().TaggerHits) {
            hTime->Fill(tHit.Time, tHit.Channel);
            hTimeToF->Fill(tHit.Time - CBTimeAvg, tHit.Channel);
            ++multiplicity[tHit.Channel];
        }
    }
    else {
        for(const auto& cand: event.Reconstructed().Candidates) {
            for(const TCluster& cluster: cand.Clusters) {
                if(cluster.DetectorType != Detector->Type)
                    continue;
                hTime->Fill(cluster.Time, cluster.CentralElement);
                ++multiplicity[cluster.CentralElement];
                const double tof = Detector->GetTimeOfFlight(cluster.Time,
                                                             cluster.CentralElement,
                                                             CBTimeAvg);
                hTimeToF->Fill(tof, cluster.CentralElement);
                for(const auto& taggerhit : event.Reconstructed().TaggerHits) {
                    const double relative_time = cluster.Time - taggerhit.Time;
                    hTimeToTagger->Fill(relative_time, cluster.CentralElement);
                }
            }
        }
    }
    for(auto& it_mult : multiplicity) {
        hTimeMultiplicity->Fill(it_mult.second, it_mult.first);
    }
}

void Time::ShowResult()
{
    canvas(GetName())
            << drawoption("colz") << hTime
            << drawoption("colz") << hTimeToTagger
            << hCBTriggerTiming
            << drawoption("colz") << hTimeToF
            << drawoption("colz") << hTimeMultiplicity
            << endc;
}

namespace ant {
namespace analysis {
namespace physics {

struct EPT_Time : Time {
    EPT_Time(const std::string& name, OptionsPtr opts) :
        Time(Detector_t::Type_t::EPT, name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(EPT_Time)

struct Tagger_Time : Time {
    Tagger_Time(const std::string& name, OptionsPtr opts) :
        Time(Detector_t::Type_t::Tagger, name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(Tagger_Time)

struct CB_Time : Time {
    CB_Time(const std::string& name, OptionsPtr opts) :
        Time(Detector_t::Type_t::CB, name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(CB_Time)

struct PID_Time : Time {
    PID_Time(const std::string& name, OptionsPtr opts) :
        Time(Detector_t::Type_t::PID, name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(PID_Time)

struct TAPS_Time : Time {
    TAPS_Time(const std::string& name, OptionsPtr opts) :
        Time(Detector_t::Type_t::TAPS, name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(TAPS_Time)

struct TAPSVeto_Time : Time {
    TAPSVeto_Time(const std::string& name, OptionsPtr opts) :
        Time(Detector_t::Type_t::TAPSVeto, name, opts)
    {}
};
AUTO_REGISTER_PHYSICS(TAPSVeto_Time)


}}}








