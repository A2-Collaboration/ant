#include "Time.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/TAPS.h"

#include "calibration/converters/MultiHitReference.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

Time::Time(const Detector_t::Type_t& detectorType,
           const string& name, OptionsPtr opts)
    :
      Physics(name, opts),
      Detector(ExpConfig::Setup::GetDetector(detectorType)),
      isTagger(dynamic_pointer_cast<TaggerDetector_t, Detector_t>(Detector) != nullptr)
{

    string detectorName(Detector_t::ToString(Detector->Type));

    // handle tagger differently
    int NrChannels = Detector->GetNChannels();
    const BinSettings TimeBins = (isTagger && NrChannels==352) ? BinSettings::RoundToBinSize(BinSettings(2000,-400,400), calibration::converter::Gains::CATCH_TDC)
                               : (isTagger) ? BinSettings(800,-400.0,400.0)
                               : BinSettings(2000,-400,400);

    hTime = HistFac.makeTH2D(detectorName + " - Time",
                             "time [ns]",
                             detectorName + " channel",
                             TimeBins,
                             BinSettings(Detector->GetNChannels()),
                             "Time"
                             );
    const AxisSettings bins_timeZoomed("t / ns", BinSettings::RoundToBinSize(BinSettings(1000,-65,65), calibration::converter::Gains::CATCH_TDC));
    hTimeToTriggerRef = HistFac.makeTH2D(
                            detectorName + " - Time relative to TriggerRef",
                            bins_timeZoomed,
                            {detectorName + " channel", Detector->GetNChannels()},
                            "hTimeToTriggerRef" // should be used for TAPS_ToF offsets...
                            );
    hTimeZoomed = HistFac.makeTH2D(
                            detectorName + " - Time (zoomed)",
                            bins_timeZoomed,
                            {detectorName + " channel", Detector->GetNChannels()},
                            "hTimeZoomed" // should be used for TAPS_ToF offsets...
                            );
    hTimeToTagger = HistFac.makeTH2D(
                        detectorName + " - Time relative to tagger",
                        "time [ns]",
                        detectorName + " channel",
                        BinSettings(2000,-1500,1500),
                        BinSettings(Detector->GetNChannels()),
                        "hTimeToTagger"
                        );
    hTriggerRefTiming = HistFac.makeTH1D("CB - Trigger timing",
                                "time [ns]",
                                "#",
                                BinSettings(500,-15,15),
                                "hTriggerRefTiming");
    hTimeMultiplicity = HistFac.makeTH2D(detectorName + " - Time Hit Multiplicity",
                                         "multiplicity",
                                         detectorName + " channel",
                                         BinSettings(8),
                                         BinSettings(Detector->GetNChannels()),
                                         "hTimeMultiplicity"
                                         );
}

void Time::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

    const double TriggerRefTime = triggersimu.GetRefTiming();
    hTriggerRefTiming->Fill(TriggerRefTime);

    // handle Tagger differently
    if(isTagger) {
	// only use 1 events to reduce background and help identify prompt peaks in early channels
	// only use clusters with more than two crystals hit (more likely to be photons)
	if (event.Reconstructed().Candidates.size() == 1 ){
		for(const auto& candidate : event.Reconstructed().Candidates){
			if (candidate.ClusterSize > 1){
        			for (const auto& tHit: event.Reconstructed().TaggerHits) {
                			hTime->Fill(tHit.Time, tHit.Channel);
                			hTimeZoomed->Fill(tHit.Time, tHit.Channel);
                			hTimeToTriggerRef->Fill(tHit.Time - TriggerRefTime, tHit.Channel);
				}
        		}
		}	
	}
    }
    else {
    	for(const auto& cand: event.Reconstructed().Candidates) {
     	   for(const TCluster& cluster: cand.Clusters) {
                if(cluster.DetectorType != Detector->Type)
                    continue;
                hTime->Fill(cluster.Time, cluster.CentralElement);
                hTimeZoomed->Fill(cluster.Time, cluster.CentralElement);
                hTimeToTriggerRef->Fill(cluster.Time - TriggerRefTime, cluster.CentralElement);
                for(const auto& taggerhit : event.Reconstructed().TaggerHits) {
                    const double relative_time = cluster.Time - taggerhit.Time;
                    hTimeToTagger->Fill(relative_time, cluster.CentralElement);
                }
            }
        }
    }

    // timing multiplicities
    for(const TDetectorReadHit& dethit : event.Reconstructed().DetectorReadHits) {
        if(dethit.DetectorType != Detector->Type)
            continue;
        if(dethit.ChannelType != Channel_t::Type_t::Timing)
            continue;
        hTimeMultiplicity->Fill(dethit.Values.size(), dethit.Channel);
    }
}

void Time::ShowResult()
{
    canvas(GetName())
            << drawoption("colz")
            << hTime
            << hTimeToTagger
            << hTriggerRefTiming
            << hTimeZoomed
            << hTimeToTriggerRef
            << hTimeMultiplicity
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








