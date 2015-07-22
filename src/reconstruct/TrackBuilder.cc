#include "TrackBuilder.h"

using namespace ant;
using namespace std;
using namespace ant::reconstruct;

void TrackBuilder::Build(std::map<Detector_t::Type_t, std::list<TCluster> >&& sorted_clusters,
                         unique_ptr<TEvent>& event)
{
    // super-stupid track matching for now
    // each track contains just one cluster

    for(const auto& it_clusters : sorted_clusters) {
        //const Detector_t::Type_t detectortype = it_clusters.first;
        const list<TCluster>& clusters = it_clusters.second;

        for(const TCluster& cluster : clusters) {
            event->Tracks.emplace_back(
                        cluster.Energy,
                        0, // time unknown...
                        cluster.Position.Theta(),
                        cluster.Position.Phi(),
                        std::vector<TCluster>{cluster}
                        );
        }
    }
}

