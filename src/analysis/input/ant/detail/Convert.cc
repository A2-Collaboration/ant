#include "Convert.h"
#include "data/Event.h"

#include "base/Detector_t.h"

#include "tree/TEvent.h"
#include "tree/TCluster.h"
#include "tree/TCandidate.h"
#include "tree/TTaggerHit.h"
#include "tree/TDataRecord.h"
#include "data/Event.h"

#include "base/Detector_t.h"

#include <memory>
#include <limits>

using namespace ant;
using namespace ant::analysis::input;
using namespace ant::analysis::data;
using namespace std;

template <class Cont1, class Cont2>
void Copy(const Cont1& from, Cont2& to) {
    to.reserve(from.size());
    for(auto& from_element : from) {
        to.emplace_back(Converter::Convert(from_element));
    }
}

Event Converter::Convert(const TEvent &event)
{
    Event antevent;

    Copy(event.Candidates,     antevent.Reconstructed().Candidates());
    Copy(event.Tagger.Hits,    antevent.Reconstructed().TaggerHits());
    Copy(event.AllClusters, antevent.Reconstructed().AllClusters());

    /// @todo The multiplicity is a much harder business, see acqu/root/src/TA2BasePhysics.cc
    /// the code there might only apply to the old trigger system before 2012
    /// so it might be best to implement such algorithms with some nicely designed interface into the
    /// pseudo-detector Trigger in expconfig/detectors

    double Esum = 0.0;
    for(const TCandidate& candidate : event.Candidates) {
        for(const TCluster& cluster: candidate.Clusters) {
            if(cluster.GetDetectorType() == Detector_t::Type_t::CB) {
                Esum += cluster.Energy;
            }
        }
    }

    auto& triggerinfos = antevent.Reconstructed().TriggerInfos();
    triggerinfos.EventID() = event.ID;
    triggerinfos.CBEenergySum() = Esum;

    return antevent;
}


shared_ptr<Candidate> Converter::Convert(const TCandidate &candidate)
{
    std::shared_ptr<Candidate> antCandidate = make_shared<Candidate>(
                        candidate.Energy,
                        candidate.Theta,
                        candidate.Phi,
                        candidate.Time,
                        0,
                        Detector_t::Any_t::None,
                        candidate.VetoEnergy,
                        candidate.TrackerEnergy
                        );

    auto det = Detector_t::Any_t::None;

    for(const TCluster& cluster: candidate.Clusters) {

        antCandidate->Clusters.emplace_back( Convert(cluster) );
        const auto& antCluster = antCandidate->Clusters.back();

        det |= antCluster.Detector;

        if(cluster.GetDetectorType() == Detector_t::Type_t::CB || cluster.GetDetectorType() == Detector_t::Type_t::TAPS) {
            antCandidate->clusterSize = cluster.Hits.size();
        }
    }

    antCandidate->detector    = det;

    return antCandidate;
}

shared_ptr<TaggerHit> Converter::Convert(const TTaggerHit& taggerhit)
{
    // careful for normal Tagger: only first electron is treated
    auto anttaggerhit = make_shared<TaggerHit>(
                            taggerhit.Electrons.front().Key,
                            taggerhit.PhotonEnergy,
                            taggerhit.Time
                            );
    return anttaggerhit;
}


Cluster Converter::Convert(const TCluster& cluster)
{

    Cluster cl(
                cluster.Energy,
                cluster.Time,
                cluster.GetDetectorType(),
                cluster.CentralElement,
                cluster.Position
                );

    if(cluster.HasFlag(TCluster::Flags_t::Split)) cl.flags.Set(Cluster::Flag::Split);
    if(cluster.HasFlag(TCluster::Flags_t::TouchesHole)) cl.flags.Set(Cluster::Flag::TouchesHole);

    for(const auto& hit : cluster.Hits) {
       Cluster::Hit anthit;
       anthit.Channel = hit.Channel;
       for(const auto& datum : hit.Data) {
           anthit.Data.emplace_back(static_cast<Channel_t::Type_t>(datum.Type), datum.Value);
       }
       cl.Hits.emplace_back(move(anthit));
    }
    return cl;
}

