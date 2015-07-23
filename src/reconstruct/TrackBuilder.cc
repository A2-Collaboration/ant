#include "TrackBuilder.h"
#include "base/Logger.h"
#include "tree/TCluster.h"
#include "expconfig/detectors/CB.h"
#include "expconfig/detectors/PID.h"
#include "expconfig/detectors/TAPS.h"

using namespace ant;
using namespace std;
using namespace ant::reconstruct;
using namespace ant::expconfig;

void TrackBuilder::Build_PID_CB(std::map<Detector_t::Type_t, std::list<TCluster> >& sorted_clusters, TEvent::tracks_t& tracks)
{
    auto& cb_clusters  = sorted_clusters[Detector_t::Type_t::CB];
    auto& pid_clusters = sorted_clusters[Detector_t::Type_t::PID];

    auto pid_cluster = pid_clusters.begin();

    while(pid_cluster != pid_clusters.end()) {

        const auto pid_phi = pid_cluster->Position.Phi();

        if(pid_cluster->Hits.size() == 1) {

            bool matched = false;

            const auto dphi = pid->dPhi(pid_cluster->Hits.at(0).Channel) / 2.0;

            auto cb_cluster = cb_clusters.begin();

            while(cb_cluster != cb_clusters.end()) {
                const auto cb_phi = cb_cluster->Position.Phi();

                if( fabs(cb_phi - pid_phi) < dphi ) { // match!

                    tracks.emplace_back(
                                cb_cluster->Energy,
                                0, // time unknown...
                                cb_cluster->Position.Theta(),
                                cb_cluster->Position.Phi(),
                                std::vector<TCluster>{*cb_cluster, *pid_cluster},
                                pid_cluster->Energy
                                );
                    cb_cluster = cb_clusters.erase(cb_cluster);
                }
            }

            if(matched) {
                pid_cluster = pid_clusters.erase(pid_cluster);
            } else {
                ++pid_cluster;
            }

        } else {
            VLOG(3) << "Strange PID Cluster encountered";
            pid_cluster = pid_clusters.erase(pid_cluster);
        }

    }
}

void TrackBuilder::Build_TAPS_Veto(std::map<Detector_t::Type_t, std::list<TCluster> >& sorted_clusters, TEvent::tracks_t& tracks)
{

}

void TrackBuilder::Catchall(std::map<Detector_t::Type_t, std::list<TCluster> >& sorted_clusters, TEvent::tracks_t& tracks)
{
    for(auto& cluster_list : sorted_clusters ) {

        const auto& detector_type = cluster_list.first;
        const auto& clusters      = cluster_list.second;

        if(detector_type == Detector_t::Type_t::PID /*|| detector_type == Detector_t::Type_t::Veto*/) {
            for(auto& c : clusters) {
                tracks.emplace_back(
                            0,
                            0, // time unknown...
                            c.Position.Theta(),
                            c.Position.Phi(),
                            std::vector<TCluster>{c},
                            c.Energy
                            );
            }
        } else if(detector_type == Detector_t::Type_t::CB || detector_type == Detector_t::Type_t::TAPS) {
            for(auto& c : clusters) {
                tracks.emplace_back(
                            c.Energy,
                            0, // time unknown...
                            c.Position.Theta(),
                            c.Position.Phi(),
                            std::vector<TCluster>{c},
                            0
                            );
            }
        } else if(detector_type == Detector_t::Type_t::MWPC0 || detector_type == Detector_t::Type_t::MWPC1 || detector_type == Detector_t::Type_t::Cherenkov) {
            for(auto& c : clusters) {
                tracks.emplace_back(
                            0,
                            0, // time unknown...
                            c.Position.Theta(),
                            c.Position.Phi(),
                            std::vector<TCluster>{c},
                            0,
                            c.Energy
                            );
            }
        }
    }
}

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

    try {
        taps = dynamic_pointer_cast<detector::TAPS>(sorted_detectors.at(Detector_t::Type_t::TAPS));
    } catch (...) {
        VLOG(3) << "Detector TAPS not initialized";
    }

//    try {
//        veto = dynamic_pointer_cast<detector::TAPSVeto>(sorted_detectors.at(Detector_t::Type_t::TAPSVeto));
//    } catch (...) {
//        VLOG(3) << "Detector TAPSVeto not initialized";
//    }

}

void TrackBuilder::Build(std::map<Detector_t::Type_t, std::list<TCluster> >&& sorted_clusters,
                         TEvent::tracks_t& tracks)
{

    if(cb && pid)
        Build_PID_CB(sorted_clusters, tracks);

    if(taps /*&& tapsveto*/)
        Build_TAPS_Veto(sorted_clusters, tracks);


    Catchall(sorted_clusters, tracks);

}

