#pragma once

#include "expconfig/ExpConfig.h"
#include "reconstruct/Reconstruct_traits.h"

#include <memory>
#include <vector>

namespace ant {

struct ClusterDetector_t;

namespace reconstruct {

class Clustering {
public:
    Clustering(const std::shared_ptr<ExpConfig::Reconstruct>& config);

    void Build(const std::shared_ptr<ClusterDetector_t>& clusterdetector,
               const TClusterHitList& clusterhits,
               TClusterList& clusters
               );
    virtual ~Clustering() = default;
};


}} // namespace ant::reconstruct
