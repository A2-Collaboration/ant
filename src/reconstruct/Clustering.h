#pragma once

#include "expconfig/ExpConfig.h"
#include "reconstruct/Reconstruct_traits.h"

#include <memory>
#include <vector>

namespace ant {

struct ClusterDetector_t;

namespace reconstruct {

class Clustering_NextGen : public Clustering_traits {
public:

    Clustering_NextGen() = default;

    virtual void Build(const ClusterDetector_t& clusterdetector,
                       const TClusterHitList& clusterhits,
                       TClusterList& clusters
                       ) const override;

    virtual ~Clustering_NextGen() = default;

};


}} // namespace ant::reconstruct
