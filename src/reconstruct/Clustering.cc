#include "Clustering.h"
#include "detail/Clustering_NextGen.h"

#include "expconfig/Detector_t.h"

#include "tree/TCluster.h"

#include "base/Logger.h"
#include "base/std_ext.h"

using namespace std;
using namespace ant;
using namespace ant::reconstruct;

HitWithEnergy_t::HitWithEnergy_t(const TDetectorReadHit* readhit,
                                              const vector<TClusterHitDatum>&& data) :
    Hit(std_ext::make_unique<TClusterHit>(readhit->Channel, data)),
    Energy(numeric_limits<double>::quiet_NaN())
{
    MaybeSetEnergy(readhit);
}

void HitWithEnergy_t::MaybeSetEnergy(const TDetectorReadHit *readhit) {
    if(readhit->GetChannelType() != Channel_t::Type_t::Integral)
        return;
    if(readhit->Values.size() != 1)
        return;
    if(isfinite(Energy)) {
        LOG(WARNING) << "Found hit in channel with more than one energy information";
        return;
    }
    Energy = readhit->Values[0];
}

void Clustering::Build(
        const shared_ptr<ClusterDetector_t>& clusterdetector,
        const list<HitWithEnergy_t>& clusterhits,
        list<TCluster>& clusters)
{
    // clustering detector, so we need additional information
    // to build the crystals_t
    list<clustering::crystal_t> crystals;
    for(const HitWithEnergy_t& clusterhit : clusterhits) {
        const auto& hit = clusterhit.Hit;
        // ignore hits without energy information
        if(!isfinite(clusterhit.Energy))
            continue;
        crystals.emplace_back(
                    clusterhit.Energy,
                    clusterdetector->GetClusterElement(hit->Channel),
                    hit.get()
                    );

    }
    // do the clustering
    vector< vector< clustering::crystal_t> > crystal_clusters;
    do_clustering(crystals, crystal_clusters);

    // now calculate some cluster properties,
    // and create TCluster out of it
    for(const auto& cluster : crystal_clusters) {
        const double cluster_energy = clustering::calc_total_energy(cluster);
        /// \todo already cut low-energy clusters here?
        TVector3 weightedPosition(0,0,0);
        double weightedSum = 0;
        std::vector<TClusterHit> clusterhits;
        clusterhits.reserve(cluster.size());
        for(const clustering::crystal_t& crystal : cluster) {
            double wgtE = clustering::calc_energy_weight(crystal.Energy, cluster_energy);
            weightedPosition += crystal.Element->Position * wgtE;
            weightedSum += wgtE;
            clusterhits.emplace_back(*crystal.Hit);
        }
        weightedPosition *= 1.0/weightedSum;
        clusters.emplace_back(
                    weightedPosition,
                    cluster_energy,
                    clusterdetector->Type,
                    clusterhits
                    );
    }
}
