#pragma once

#include <memory>
#include <list>
#include <vector>

namespace ant {

class ClusterDetector_t;
class TDetectorReadHit;
class TCluster;
class TClusterHit;
class TClusterHitDatum;

namespace reconstruct {

/**
 * @brief The HitWithEnergy_t struct
 *
 * Stores a TClusterHit together with its energy, if
 * exactly one TDetectorReadHit was of channel type integral
 *
 */
struct HitWithEnergy_t {
    std::unique_ptr<TClusterHit> Hit;
    double Energy;
    void MaybeSetEnergy(const TDetectorReadHit* readhit);
    HitWithEnergy_t(const TDetectorReadHit* readhit,
                    const std::vector<TClusterHitDatum>&& data);
};

class Clustering {
public:
    void Build(const std::shared_ptr<ClusterDetector_t>& clusterdetector,
               const std::list<HitWithEnergy_t>& clusterhits,
               std::list<TCluster>& clusters
               );
};


}} // namespace ant::reconstruct
