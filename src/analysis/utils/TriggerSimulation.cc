#include "TriggerSimulation.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;

bool TriggerSimulation::ProcessEvent(const TEvent& event)
{
    info.Reset();

    const auto& recon = event.Reconstructed();

    const bool isMC = recon.ID.isSet(TID::Flags_t::MC);

    // CBTiming as average over cluster timings
    {
        if(isMC) {
            info.CBTiming = 0;
        }
        else {
            double TimeEsum = 0.0;
            double TimeE = 0.0;
            for(const auto& cluster : recon.Clusters) {
                if(cluster.DetectorType != Detector_t::Type_t::CB)
                    continue;
                // ignore weird clusters
                if(!cluster.isSane())
                    continue;
                TimeEsum += cluster.Energy;
                TimeE += cluster.Energy*cluster.Time;
            }
            info.CBTiming = TimeE/TimeEsum;
        }
    }

    // CBEsum is sum over detector read hits
    {
        info.CBEnergySum = 0.0;
        for(const TDetectorReadHit& dethit : recon.DetectorReadHits) {
            if(dethit.DetectorType != Detector_t::Type_t::CB)
                continue;
            if(dethit.ChannelType != Channel_t::Type_t::Integral)
                continue;
            // one could also use uncalibrated values
            // with some fixed constant? or average from all gains?
            for(auto& energy : dethit.Values)
                info.CBEnergySum += energy.Calibrated;
        }
    }

    /// \todo Improve the trigger simulation on MC here, could also depend on config given to
    /// TriggerSimulation ctor
    // assume data has always triggered, no simulation needed
    info.hasTriggered = !isMC || info.CBEnergySum > 550;

    /// \todo The multiplicity is a much harder business, see acqu/root/src/TA2BasePhysics.cc
    /// the code there might only apply to the old trigger system before 2012

    // return true if information is complete and sane
    return info.IsSane();
}
