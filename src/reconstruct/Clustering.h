#pragma once

#include <memory>
#include <list>
#include <vector>

#include "expconfig/ExpConfig.h"

namespace ant {

class ClusterDetector_t;
class TDetectorReadHit;
class TCluster;
class TClusterHit;
class TClusterHitDatum;

namespace reconstruct {

/**
 * @brief The AdaptorTClusterHit struct
 *
 * Stores a TClusterHit together with its energy and timing
 *
 */
struct AdaptorTClusterHit {
    std::unique_ptr<TClusterHit> Hit;
    double Energy;
    double Time;
    void SetFields(const TDetectorReadHit* readhit);
    AdaptorTClusterHit(const TDetectorReadHit* readhit,
                       const std::vector<TClusterHitDatum>&& data);
};

class Clustering {
public:
    Clustering(const std::shared_ptr<ExpConfig::Reconstruct>& config);

    void Build(const std::shared_ptr<ClusterDetector_t>& clusterdetector,
               const std::list<AdaptorTClusterHit>& clusterhits,
               std::list<TCluster>& clusters
               );
    virtual ~Clustering() = default;
};


}} // namespace ant::reconstruct
