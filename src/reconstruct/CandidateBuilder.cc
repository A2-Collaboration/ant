#include "CandidateBuilder.h"
#include "base/Logger.h"
#include "tree/TCluster.h"
#include "expconfig/detectors/CB.h"
#include "expconfig/detectors/PID.h"
#include "expconfig/detectors/TAPS.h"
#include "expconfig/detectors/TAPSVeto.h"
#include "base/std_ext/math.h"

#include "TVector2.h"


using namespace ant;
using namespace std;
using namespace ant::reconstruct;
using namespace ant::expconfig;

CandidateBuilder::CandidateBuilder(
        const CandidateBuilder::sorted_detectors_t& sorted_detectors,
        const std::shared_ptr<ExpConfig::Reconstruct>& _config
        ) :
    config(_config->GetCandidateBuilderConfig())
{

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

void CandidateBuilder::Build_PID_CB(std::map<Detector_t::Type_t, std::list<TCluster> >& sorted_clusters,
                                    TEvent::candidates_t& candidates,
                                    std::vector<TCluster>& all_clusters)
{
    auto it_cb_clusters = sorted_clusters.find(Detector_t::Type_t::CB);
    if(it_cb_clusters == sorted_clusters.end())
        return;
    auto& cb_clusters  = it_cb_clusters->second;
    if(cb_clusters.empty())
        return;

    auto it_pid_clusters = sorted_clusters.find(Detector_t::Type_t::PID);
    if(it_pid_clusters == sorted_clusters.end())
        return;
    auto& pid_clusters  = it_pid_clusters->second;
    if(pid_clusters.empty())
        return;

    auto pid_cluster = pid_clusters.begin();

    const auto dphi_max = (pid->dPhi(pid_cluster->Hits.at(0).Channel) + config.PID_Phi_Epsilon) / 2.0;

    while(pid_cluster != pid_clusters.end()) {

        const auto pid_phi = pid_cluster->Position.Phi();

        if(pid_cluster->Hits.size() == 1) {

            bool matched = false;

            auto cb_cluster = cb_clusters.begin();

            while(cb_cluster != cb_clusters.end()) {
                const auto cb_phi = cb_cluster->Position.Phi();

                // calculate phi angle difference.
                // Phi_mpi_pi() takes care of wrap-arounds at 180/-180 deg
                const auto dphi = fabs(TVector2::Phi_mpi_pi(cb_phi - pid_phi));
                if(dphi < dphi_max ) { // match!

                    candidates.emplace_back(
                                cb_cluster->Energy,
                                cb_cluster->Time,
                                cb_cluster->Position.Theta(),
                                cb_cluster->Position.Phi(),
                                std::vector<TCluster>{*cb_cluster, *pid_cluster},
                                pid_cluster->Energy
                                );
                    all_clusters.emplace_back(move(*cb_cluster));
                    cb_cluster = cb_clusters.erase(cb_cluster);
                    matched = true;
                }
                else {
                    ++cb_cluster;
                }
            }

            if(matched) {
                all_clusters.emplace_back(move(*pid_cluster));
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

void CandidateBuilder::Build_TAPS_Veto(std::map<Detector_t::Type_t, std::list<TCluster> >& sorted_clusters,
                                       TEvent::candidates_t& candidates,
                                       std::vector<TCluster>& all_clusters)
{
    auto it_taps_clusters = sorted_clusters.find(Detector_t::Type_t::TAPS);
    if(it_taps_clusters == sorted_clusters.end())
        return;
    auto& taps_clusters  = it_taps_clusters->second;
    if(taps_clusters.empty())
        return;

    auto it_veto_clusters = sorted_clusters.find(Detector_t::Type_t::TAPSVeto);
    if(it_veto_clusters == sorted_clusters.end())
        return;
    auto& veto_clusters  = it_veto_clusters->second;
    if(veto_clusters.empty())
        return;

    const auto element_radius2 = std_ext::sqr(tapsveto->GetElementRadius());

    auto veto_cluster = veto_clusters.begin();
    while(veto_cluster != veto_clusters.end()) {

        bool matched = false;

        const TVector3& vpos = veto_cluster->Position;

        auto taps_cluster = taps_clusters.begin();

        while(taps_cluster != taps_clusters.end()) {

            const TVector3& tpos = taps_cluster->Position;
            const TVector3 d = tpos - vpos;

            if( d.XYvector().Mod() < element_radius2 ) {
                candidates.emplace_back(
                            taps_cluster->Energy,
                            taps_cluster->Time,
                            taps_cluster->Position.Theta(),
                            taps_cluster->Position.Phi(),
                            std::vector<TCluster>{*taps_cluster, *veto_cluster},
                            veto_cluster->Energy
                            );
                all_clusters.emplace_back(move(*taps_cluster));
                taps_cluster = taps_clusters.erase(taps_cluster);
                matched = true;
            } else {
                ++taps_cluster;
            }
        }

        if(matched) {
            all_clusters.emplace_back(move(*veto_cluster));
            veto_cluster = veto_clusters.erase(veto_cluster);
        } else {
            ++veto_cluster;
        }
    }
}

void CandidateBuilder::Catchall(std::map<Detector_t::Type_t, std::list<TCluster> >& sorted_clusters,
                                TEvent::candidates_t& candidates,
                                std::vector<TCluster>& all_clusters)
{
    for(auto& cluster_list : sorted_clusters ) {

        const auto& detector_type = cluster_list.first;
        auto& clusters      = cluster_list.second;

        if(option_allowSingleVetoClusters &&
           (detector_type == Detector_t::Type_t::PID || detector_type == Detector_t::Type_t::TAPSVeto)) {
            for(auto& c : clusters) {
                candidates.emplace_back(
                            0,
                            c.Time,
                            c.Position.Theta(),
                            c.Position.Phi(),
                            std::vector<TCluster>{c},
                            c.Energy
                            );
                all_clusters.emplace_back(move(c));
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
                all_clusters.emplace_back(move(c));
            }
            clusters.clear();
        } else if(detector_type == Detector_t::Type_t::MWPC0 || detector_type == Detector_t::Type_t::MWPC1
                  || detector_type == Detector_t::Type_t::Cherenkov) {
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
                all_clusters.emplace_back(move(c));
            }
            clusters.clear();
        }
    }
}





void CandidateBuilder::BuildCandidates(
        std::map<Detector_t::Type_t, std::list<TCluster> >& sorted_clusters,
        TEvent::candidates_t& candidates,
        std::vector<TCluster>& all_clusters
        )
{
    // build candidates
    if(cb && pid)
        Build_PID_CB(sorted_clusters, candidates, all_clusters);

    if(taps && tapsveto)
        Build_TAPS_Veto(sorted_clusters, candidates, all_clusters);

    Catchall(sorted_clusters, candidates, all_clusters);
}

void CandidateBuilder::Build(
        std::map<Detector_t::Type_t, std::list<TCluster> > sorted_clusters,
        TEvent::candidates_t& candidates,
        std::vector<TCluster>& all_clusters
        )
{
    // search for clusters which are not sane or don't pass thresholds
    // add them to all_clusters with unmatched flag set
    for(auto& det_entry : sorted_clusters) {
        auto detectortype = det_entry.first;
        auto& clusters = det_entry.second;
        auto it_cluster = clusters.begin();
        while(it_cluster != clusters.end()) {
            double threshold = 0.0;
            if(detectortype == Detector_t::Type_t::CB)
                threshold = config.CB_ClusterThreshold;
            if(detectortype == Detector_t::Type_t::TAPS)
                threshold = config.TAPS_ClusterThreshold;

            // do not remove clusters which are sane and pass the thresholds
            if(it_cluster->isSane() && it_cluster->Energy > threshold) {
                ++it_cluster;
                continue;
            }


            it_cluster->SetFlag(TCluster::Flags_t::Unmatched);
            all_clusters.emplace_back(move(*it_cluster));
            it_cluster = clusters.erase(it_cluster);
        }
    }

    BuildCandidates(sorted_clusters, candidates, all_clusters);

    // add remaining unmatched clusters to all_clusters with Unmatched flag set
    for(auto& det_entry : sorted_clusters) {
        for(auto& cluster : det_entry.second) {
            cluster.SetFlag(TCluster::Flags_t::Unmatched);
            all_clusters.emplace_back(move(cluster));
        }
    }
}

