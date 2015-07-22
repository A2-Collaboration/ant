#include "TrackBuilder.h"
#include "base/Logger.h"
#include "tree/TCluster.h"
#include "expconfig/detectors/CB.h"
#include "expconfig/detectors/PID.h"

using namespace ant;
using namespace std;
using namespace ant::reconstruct;
using namespace ant::expconfig;

TrackBuilder::TrackBuilder(const TrackBuilder::sorted_detectors_t& sorted_detectors) {

    try {
        cb  = dynamic_pointer_cast<detector::CB>(sorted_detectors.at(Detector_t::Type_t::CB));
    } catch (...) {
        VLOG(3) << "Detector CB not initialized";
    }

    try {
        pid = dynamic_pointer_cast<detector::PID>(sorted_detectors.at(Detector_t::Type_t::PID));
    } catch (...) {
        VLOG(3) << "Detector PID not initialized";
    }

}

void TrackBuilder::Build(std::map<Detector_t::Type_t, std::list<TCluster> >&& sorted_clusters,
                         unique_ptr<TEvent>& event)
{

    if(cb) {

        for(TCluster& cb_cluster : sorted_clusters[Detector_t::Type_t::CB]) {
            const auto cb_phi = cb_cluster.Position.Phi();

            event->Tracks.emplace_back(
                        cb_cluster.Energy,
                        0, // time unknown...
                        cb_cluster.Position.Theta(),
                        cb_cluster.Position.Phi(),
                        std::vector<TCluster>{cb_cluster}
                        );

            TTrack& track = event->Tracks.back();

            if(pid) {
                for(TCluster& pid_cluster : sorted_clusters[Detector_t::Type_t::PID]) {
                    const auto pid_phi = pid_cluster.Position.Phi();

                    if(pid_cluster.Hits.size() == 1) {

                        const auto dphi = pid->dPhi(pid_cluster.Hits.at(0).Channel) / 2.0;

                        if( fabs(cb_phi - pid_phi) < dphi ) { // match!
                            track.Clusters.emplace_back(pid_cluster);
                            track.VetoEnergy = pid_cluster.Energy;
                        }

                    } else {
                        VLOG(3) << "Empty PID Cluster encountered";
                    }

                }
            }
        }
    } /* if cb */

}

