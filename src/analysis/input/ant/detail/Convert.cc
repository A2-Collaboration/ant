#include "Convert.h"
#include <memory>
#include "data/Event.h"
#include "tree/TEvent.h"
#include "tree/TCluster.h"
#include "tree/TCandidate.h"
#include "tree/TTaggerHit.h"
#include "tree/TDataRecord.h"
#include "analysis/data/Candidate.h"
#include "analysis/data/TaggerHit.h"
#include "analysis/Detector.h"
#include "expconfig/Detector_t.h"
#include <limits>

using namespace ant;
using namespace ant::input;
using namespace std;

template <class Cont1, class Cont2>
void Copy(const Cont1& from, Cont2& to) {
    to.reserve(from.size());
    for(auto& from_element : from) {
        to.emplace_back(ant::input::Convert(from_element));
    }
}

shared_ptr<ant::Event> input::Convert(const TEvent &event)
{
    shared_ptr<Event> antevent = make_shared<Event>();

    Copy(event.Candidates,     antevent->Reconstructed().Candidates());
    Copy(event.Tagger.Hits, antevent->Reconstructed().TaggerHits());

    return std::move(antevent);
}


std::shared_ptr<Candidate> input::Convert(const TCandidate &candidate)
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

std::shared_ptr<TaggerHit> input::Convert(const TTaggerHit& taggerhit)
{
    ///@todo implement something for tagger channel
    auto anttaggerhit = make_shared<TaggerHit>(
                            0,
                            taggerhit.PhotonEnergy,
                            taggerhit.Time
                            );
    return anttaggerhit;
}


Cluster input::Convert(const TCluster& cluster)
{

    return Cluster(
                cluster.Energy,
                cluster.Time,
                detector_t(Detector_t::ToBitfield(cluster.GetDetectorType())),
                cluster.CentralElement
                );
}

