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
    Clustering();

    virtual void Build(const ClusterDetector_t& clusterdetector,
                       const TClusterHitList& clusterhits,
                       TClusterList& clusters
                       );
    virtual ~Clustering() = default;
};


}} // namespace ant::reconstruct
