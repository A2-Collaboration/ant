#include "Clustering.h"
#include "detail/Clustering_NextGen.h"

#include "base/Detector_t.h"

#include "tree/TCluster.h"
#include "tree/TDetectorReadHit.h"

#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::reconstruct;

Clustering::Clustering(const shared_ptr<ExpConfig::Reconstruct>&)
{
}

void Clustering::Build(const shared_ptr<ClusterDetector_t>& clusterdetector,
        const TClusterHitList& clusterhits,
        TClusterList& clusters)
{
    // clustering detector, so we need additional information
    // to build the crystals_t
    list<clustering::crystal_t> crystals;
    for(const TClusterHit& hit : clusterhits) {
        // ignore hits without energy or time information
        if(!hit.IsSane()) {
            continue;
        }
        crystals.emplace_back(
                    hit.Energy,
                    clusterdetector->GetClusterElement(hit.Channel),
                    addressof(hit)
                    );
    }

    // do the clustering (calls detail/Clustering_NextGen.h code)
    vector< clustering::cluster_t > crystal_clusters;
    clustering::do_clustering(crystals, crystal_clusters);

    // now calculate some cluster properties,
    // and create TCluster out of it

    for(const clustering::cluster_t& cluster : crystal_clusters) {
        const double cluster_energy = clustering::calc_total_energy(cluster);

        vec3   weightedPosition(0,0,0);
        double weightedSum = 0;

        double   cluster_time        = numeric_limits<double>::quiet_NaN();
        double   cluster_maxenergy   = 0;
        unsigned cluster_max_channel = 0;
        double   cluster_shortenergy = numeric_limits<double>::quiet_NaN();
        bool     cluster_toucheshole = false;

        std::vector<TClusterHit> clusterhits;
        clusterhits.reserve(cluster.size());

        for(const clustering::crystal_t& crystal : cluster) {
            double wgtE = clustering::calc_energy_weight(crystal.Energy, cluster_energy);
            weightedPosition += crystal.Element->Position * wgtE;
            weightedSum += wgtE;
            clusterhits.emplace_back(*crystal.Hit);
            if(cluster_maxenergy<=crystal.Energy) {
                cluster_toucheshole = crystal.Element->TouchesHole;
                cluster_time = crystal.Hit->Time;
                cluster_maxenergy = crystal.Energy;
                cluster_max_channel = crystal.Element->Channel;
                // search for short energy
                for(const TClusterHit::Datum& datum : crystal.Hit->Data) {
                    if(datum.Type == Channel_t::Type_t::IntegralShort) {
                        cluster_shortenergy = datum.Value;
                        break;
                    }
                }
            }
        }
        weightedPosition *= 1.0/weightedSum;


        clusters.emplace_back(
                    weightedPosition,
                    cluster_energy,
                    cluster_time,
                    clusterdetector->Type,
                    cluster_max_channel,
                    clusterhits
                    );
        auto& the_cluster = clusters.back();
        if(isfinite(cluster_shortenergy))
            the_cluster.ShortEnergy = cluster_shortenergy;
        if(cluster.Split)
            the_cluster.SetFlag(TCluster::Flags_t::Split);
        if(cluster_toucheshole)
            the_cluster.SetFlag(TCluster::Flags_t::TouchesHole);
    }
}

