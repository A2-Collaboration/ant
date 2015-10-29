#include "Clustering.h"
#include "detail/Clustering_NextGen.h"

#include "base/Detector_t.h"

#include "tree/TCluster.h"
#include "tree/TDetectorRead.h"

#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::reconstruct;

AdaptorTClusterHit::AdaptorTClusterHit(const TDetectorReadHit* readhit,
                                       const vector<TClusterHitDatum>&& data) :
    Hit(std_ext::make_unique<TClusterHit>(readhit->Channel, data)),
    Energy(numeric_limits<double>::quiet_NaN()), // use nan as unset indicator
    Time(numeric_limits<double>::quiet_NaN())
{
    SetFields(readhit);
}

void AdaptorTClusterHit::SetFields(const TDetectorReadHit *readhit) {
    if(std::isnan(Energy) && readhit->GetChannelType() == Channel_t::Type_t::Integral)
        Energy = readhit->Values[0];
    if(std::isnan(Time) && readhit->GetChannelType() == Channel_t::Type_t::Timing)
        Time = readhit->Values[0];
}

Clustering::Clustering(const shared_ptr<ExpConfig::Reconstruct>&)
{
}

void Clustering::Build(
        const shared_ptr<ClusterDetector_t>& clusterdetector,
        const list<AdaptorTClusterHit>& clusterhits,
        list<TCluster>& clusters)
{
    // clustering detector, so we need additional information
    // to build the crystals_t
    list<clustering::crystal_t> crystals;
    for(const AdaptorTClusterHit& clusterhit : clusterhits) {
        const auto& hit = clusterhit.Hit;
        // ignore hits without energy or time information
        if(!isfinite(clusterhit.Energy) || !isfinite(clusterhit.Time)) {
            // we're not allowed to throw anything away
            // so add some strange single hit cluster with no energy information here
            clusters.emplace_back(
                        clusterdetector->GetPosition(hit->Channel),
                        clusterhit.Energy,
                        clusterhit.Time,
                        clusterdetector->Type,
                        hit->Channel,
                        vector<TClusterHit>{*hit}
                        );
            continue;
        }
        crystals.emplace_back(
                    clusterhit.Energy,
                    clusterdetector->GetClusterElement(hit->Channel),
                    addressof(clusterhit)
                    );
    }

    // do the clustering (calls detail/Clustering_NextGen.h code)
    vector< clustering::cluster_t > crystal_clusters;
    clustering::do_clustering(crystals, crystal_clusters);

    // now calculate some cluster properties,
    // and create TCluster out of it

    for(const clustering::cluster_t& cluster : crystal_clusters) {
        const double cluster_energy = clustering::calc_total_energy(cluster);

        TVector3 weightedPosition(0,0,0);
        double weightedSum = 0;

        double   cluster_time        = numeric_limits<double>::quiet_NaN();
        double   cluster_maxenergy   = 0;
        unsigned cluster_max_channel = 0;

        std::vector<TClusterHit> clusterhits;
        clusterhits.reserve(cluster.size());

        for(const clustering::crystal_t& crystal : cluster) {
            double wgtE = clustering::calc_energy_weight(crystal.Energy, cluster_energy);
            weightedPosition += crystal.Element->Position * wgtE;
            weightedSum += wgtE;
            clusterhits.emplace_back(*(crystal.Hit->Hit));
            if(cluster_maxenergy<=crystal.Energy) {
                cluster_time = crystal.Hit->Time;
                cluster_maxenergy = crystal.Energy;
                cluster_max_channel = crystal.Element->Channel;
            }
        }
        weightedPosition *= 1.0/weightedSum;
        TCluster the_cluster(
                    weightedPosition,
                    cluster_energy,
                    cluster_time,
                    clusterdetector->Type,
                    cluster_max_channel,
                    clusterhits
                    );
        if(cluster.Split)
            the_cluster.SetFlag(TCluster::Flags_t::Split);
        clusters.emplace_back(move(the_cluster));
    }
}

