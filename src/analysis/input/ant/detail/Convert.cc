#include "Convert.h"
#include "data/Event.h"

#include "Detector.h"

#include "tree/TEvent.h"
#include "tree/TCluster.h"
#include "tree/TCandidate.h"
#include "tree/TTaggerHit.h"
#include "tree/TDataRecord.h"


#include "expconfig/Detector_t.h"

#include <memory>
#include <limits>

using namespace ant;
using namespace ant::input;
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
    Copy(event.InsaneClusters, antevent.Reconstructed().InsaneClusters());

    return antevent;
}


shared_ptr<Candidate> Converter::Convert(const TCandidate &candidate)
{
    /// @todo implement cluster size
    /// @todo add clusters to ant::Candidate
    std::shared_ptr<Candidate> antCandidate = make_shared<Candidate>(
                        candidate.Energy,
                        candidate.Theta,
                        candidate.Phi,
                        candidate.Time,
                        0,
                        ant::detector_t::None,
                        candidate.VetoEnergy,
                        candidate.TrackerEnergy
                        );

    auto det = ant::detector_t::None;

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
    ///@todo implement something for tagger channel
    auto anttaggerhit = make_shared<TaggerHit>(
                            0,
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
                detector_t(Detector_t::ToBitfield(cluster.GetDetectorType())),
                cluster.CentralElement
                );
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

