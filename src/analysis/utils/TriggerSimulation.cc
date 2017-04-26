#include "TriggerSimulation.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;

TriggerSimulation::TriggerSimulation() :
    config(ExpConfig::Setup::Get().GetTriggerSimuConfig()),
    random_CBESum_threshold(config.CBESum_Edge, config.CBESum_Width)
{}

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

    // CBEsum is sum over detector read hits (except possibly ignored channels)
    {
        info.CBEnergySum = 0.0;
        for(const TDetectorReadHit& dethit : recon.DetectorReadHits) {
            if(dethit.DetectorType != Detector_t::Type_t::CB)
                continue;
            if(dethit.ChannelType != Channel_t::Type_t::Integral)
                continue;
            if(std_ext::contains(config.CBESum_MissingElements, dethit.Channel))
                continue;
            // one could also use uncalibrated values
            // with some fixed constant? or average from all gains?
            for(auto& energy : dethit.Values)
                info.CBEnergySum += energy.Calibrated;
        }
    }

    if(isMC) {
        if(config.Type == config_t::Type_t::CBESum) {
            // lazy init random generator with timestamp of (first) event
            if(!random_gen) {
                random_gen = std_ext::make_unique<std::default_random_engine>(
                                 event.Reconstructed().ID.Timestamp
                                 );
            }
            info.hasTriggered = info.CBEnergySum > random_CBESum_threshold(*random_gen);
        }
        // may implement other trigger simulations on MC here
        else {
            info.hasTriggered = true;
            LOG_N_TIMES(1, WARNING) << "Cannot simulate trigger on MC, no config provided by setup";
        }
    }
    else {
        // assume data has always triggered, no simulation needed
        info.hasTriggered = true;
    }

    /// \todo The multiplicity is a much harder business, see acqu/root/src/TA2BasePhysics.cc
    /// the code there might only apply to the old trigger system before 2012

    // return true if information is complete and sane
    return info.IsSane();
}

double TriggerSimulation::GetCorrectedTaggerTime(const TTaggerHit& taggerhit) const {
    return taggerhit.Time - GetRefTiming();
}
