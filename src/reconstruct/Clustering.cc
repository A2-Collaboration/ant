#include "Clustering.h"
#include "detail/Clustering_NextGen.h"

#include "base/Detector_t.h"

#include "tree/TCluster.h"

using namespace std;
using namespace ant;
using namespace ant::reconstruct;

bool check_TClusterHit(const TClusterHit& hit, const ClusterDetector_t& clusterdetector) {
    if(hit.IsSane())
        return true;
    if(clusterdetector.HasElementFlags(hit.Channel, Detector_t::ElementFlag_t::BadTDC)
       && isfinite(hit.Energy))
        return true;
    return false;
}

void Clustering_NextGen::Build(const ClusterDetector_t& clusterdetector,
        const TClusterHitList& clusterhits,
        TClusterList& clusters) const
{
    // clustering detector, so we need additional information
    // to build the crystals_t
    list<clustering::crystal_t> crystals;
    for(const TClusterHit& hit : clusterhits) {
        // try to include as many hits as possible
        if(!check_TClusterHit(hit, clusterdetector)) {
            continue;
        }
        crystals.emplace_back(
                    hit.Energy,
                    clusterdetector.GetClusterElement(hit.Channel),
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

        clusters.emplace_back(
                    vec3(0,0,0),
                    cluster_energy,
                    std_ext::NaN,
                    clusterdetector.Type,
                    0
                    );
        auto& the_cluster = clusters.back();

        auto& clusterhits = the_cluster.Hits;
        clusterhits.reserve(cluster.size());

        double weightedSum = 0;
        double cluster_maxenergy = 0;
        double cluster_maxenergy_goodTDC = 0;
        double centralcrystal_time = std_ext::NaN;
        bool crystalTouchesHole = false;
        for(const clustering::crystal_t& crystal : cluster) {

            clusterhits.emplace_back(*crystal.Hit);

            double wgtE = clustering::calc_energy_weight(crystal.Energy, cluster_energy);
            the_cluster.Position += crystal.Element->Position * wgtE;
            weightedSum += wgtE;

            crystalTouchesHole |= crystal.Element->TouchesHole;

            // search for crystal with maximum energy
            // which is defined as the central element
            if(crystal.Energy >= cluster_maxenergy) {
                cluster_maxenergy = crystal.Energy;

                the_cluster.SetFlag(TCluster::Flags_t::TouchesHoleCentral, crystal.Element->TouchesHole);

                centralcrystal_time = crystal.Hit->Time;
                the_cluster.CentralElement = crystal.Element->Channel;
                // search for short energy
                for(const TClusterHit::Datum& datum : crystal.Hit->Data) {
                    if(datum.Type == Channel_t::Type_t::IntegralShort) {
                        the_cluster.ShortEnergy = datum.Value.Calibrated;
                        break;
                    }
                }
            }
            // search for the crystal with maximum energy which is NOT labeled with bad TDC
            if(crystal.Energy >= cluster_maxenergy_goodTDC && (crystal.Element->Flags != Detector_t::ElementFlag_t::BadTDC)) {
                the_cluster.Time = crystal.Hit->Time;
                cluster_maxenergy_goodTDC = crystal.Energy;
            }

        }

        the_cluster.Position *= 1.0/weightedSum;

        if(cluster.Split)
            the_cluster.SetFlag(TCluster::Flags_t::Split);

        if(crystalTouchesHole)
            the_cluster.SetFlag(TCluster::Flags_t::TouchesHoleCrystal);

        // if the cluster still has no time assigned it means it only contains BadTDC crystals
        //  then just use the time for the central element
        if(!(std::isfinite(the_cluster.Time)))
            the_cluster.Time = centralcrystal_time;
    }
}

