#include "CandidateBuilder.h"
#include "base/Logger.h"
#include "tree/TCluster.h"
#include "expconfig/detectors/CB.h"
#include "expconfig/detectors/PID.h"
#include "expconfig/detectors/TAPS.h"
#include "expconfig/detectors/TAPSVeto.h"
#include "base/std_ext.h"

#include "TVector2.h"


using namespace ant;
using namespace std;
using namespace ant::reconstruct;
using namespace ant::expconfig;

void CandidateBuilder::Build_PID_CB(std::map<Detector_t::Type_t, std::list<TCluster> >& sorted_clusters, TEvent::candidates_t& candidates)
{
    auto& cb_clusters  = sorted_clusters[Detector_t::Type_t::CB];
    auto& pid_clusters = sorted_clusters[Detector_t::Type_t::PID];

    auto pid_cluster = pid_clusters.begin();

    while(pid_cluster != pid_clusters.end()) {

        const auto pid_phi = pid_cluster->Position.Phi();

        if(pid_cluster->Hits.size() == 1) {

            bool matched = false;

            const auto dphi_max = pid->dPhi(pid_cluster->Hits.at(0).Channel) / 2.0;

            auto cb_cluster = cb_clusters.begin();

            while(cb_cluster != cb_clusters.end()) {
                const auto cb_phi = cb_cluster->Position.Phi();

                // calculate phi angle difference.
                // Phi_mpi_pi() takes care of wrap-arounds at 180/-180 deg
                const auto dphi = fabs(TVector2::Phi_mpi_pi(cb_phi - pid_phi));
                if(  dphi < dphi_max ) { // match!

                    candidates.emplace_back(
                                cb_cluster->Energy,
                                cb_cluster->Time,
                                cb_cluster->Position.Theta(),
                                cb_cluster->Position.Phi(),
                                std::vector<TCluster>{*cb_cluster, *pid_cluster},
                                pid_cluster->Energy
                                );
                    cb_cluster = cb_clusters.erase(cb_cluster);
                }
                else {
                    ++cb_cluster;
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

void CandidateBuilder::Build_TAPS_Veto(std::map<Detector_t::Type_t, std::list<TCluster> >& sorted_clusters, TEvent::candidates_t& candidates)
{
    auto& veto_clusters = sorted_clusters[Detector_t::Type_t::TAPSVeto];
    auto& taps_clusters = sorted_clusters[Detector_t::Type_t::TAPS];

    const auto element_radius2 = std_ext::sqr(tapsveto->GetElementRadius());

    auto veto = veto_clusters.begin();
    while(veto != veto_clusters.end()) {

        bool matched = false;

        const TVector3& vpos = veto->Position;

        auto taps = taps_clusters.begin();

        while(taps != taps_clusters.end()) {

            const TVector3& tpos = taps->Position;
            const TVector3 d = tpos - vpos;

            if( d.XYvector().Mod() < element_radius2 ) {
                candidates.emplace_back(
                            taps->Energy,
                            taps->Time,
                            taps->Position.Theta(),
                            taps->Position.Phi(),
                            std::vector<TCluster>{*taps, *veto},
                            veto->Energy
                            );
                taps = taps_clusters.erase(taps);
            } else {
                ++taps;
            }
        }

        if(matched) {
            veto = veto_clusters.erase(veto);
        } else {
            ++veto;
        }
    }
}

void CandidateBuilder::Catchall(std::map<Detector_t::Type_t, std::list<TCluster> >& sorted_clusters, TEvent::candidates_t& candidates)
{
    for(auto& cluster_list : sorted_clusters ) {

        const auto& detector_type = cluster_list.first;
        auto& clusters      = cluster_list.second;

        if(option_allowSingleVetoClusters && (detector_type == Detector_t::Type_t::PID || detector_type == Detector_t::Type_t::TAPSVeto)) {
            for(auto& c : clusters) {
                candidates.emplace_back(
                            0,
                            c.Time,
                            c.Position.Theta(),
                            c.Position.Phi(),
                            std::vector<TCluster>{c},
                            c.Energy
                            );
            }
            clusters.clear();
        } else if(detector_type == Detector_t::Type_t::CB || detector_type == Detector_t::Type_t::TAPS) {
            for(auto& c : clusters) {
                candidates.emplace_back(
                            c.Energy,
                            c.Time,
                            c.Position.Theta(),
                            c.Position.Phi(),
                            std::vector<TCluster>{c},
                            0
                            );
            }
            clusters.clear();
        } else if(detector_type == Detector_t::Type_t::MWPC0 || detector_type == Detector_t::Type_t::MWPC1 || detector_type == Detector_t::Type_t::Cherenkov) {
            for(auto& c : clusters) {
                candidates.emplace_back(
                            0,
                            c.Time,
                            c.Position.Theta(),
                            c.Position.Phi(),
                            std::vector<TCluster>{c},
                            0,
                            c.Energy
                            );
            }
            clusters.clear();
        }
    }
}

CandidateBuilder::CandidateBuilder(const CandidateBuilder::sorted_detectors_t& sorted_detectors) {

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

    try {
        tapsveto = dynamic_pointer_cast<detector::TAPSVeto>(sorted_detectors.at(Detector_t::Type_t::TAPSVeto));
    } catch (...) {
        VLOG(3) << "Detector TAPSVeto not initialized";
    }

}

void CandidateBuilder::Build(std::map<Detector_t::Type_t, std::list<TCluster> >&& sorted_clusters,
                         TEvent::candidates_t& tracks, std::vector<TCluster>& insane_clusters)
{

    if(cb && pid)
        Build_PID_CB(sorted_clusters, tracks);

    if(taps && tapsveto)
        Build_TAPS_Veto(sorted_clusters, tracks);

    Catchall(sorted_clusters, tracks);

    // move the rest to insane clusters
    for(auto& det_entry : sorted_clusters) {
        for(auto& cluster : det_entry.second) {
            insane_clusters.emplace_back(cluster);
        }
    }
}

