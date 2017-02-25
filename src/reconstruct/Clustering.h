#pragma once

#include "expconfig/ExpConfig.h"
#include "reconstruct/Reconstruct_traits.h"

#include <memory>
#include <vector>

namespace ant {

struct ClusterDetector_t;

namespace reconstruct {

class Clustering_traits {
public:
    virtual void Build(const ClusterDetector_t& clusterdetector,
                       const TClusterHitList& clusterhits,
                       TClusterList& clusters
                       ) const = 0;
    virtual ~Clustering_traits() = default;
};

class Clustering_NextGen : public Clustering_traits {
public:

    Clustering_NextGen() = default;

    virtual void Build(const ClusterDetector_t& clusterdetector,
                       const TClusterHitList& clusterhits,
                       TClusterList& clusters
                       ) const override;

    virtual ~Clustering_NextGen() = default;

};

class Clustering_Sergey : public Clustering_traits {
public:

    Clustering_Sergey();

    virtual void Build(const ClusterDetector_t& clusterdetector,
                       const TClusterHitList& clusterhits,
                       TClusterList& clusters
                       ) const override;

    virtual ~Clustering_Sergey();

private:
    // PIMPL idiom
    struct Impl;
    const std::unique_ptr<Impl> impl;
};


}} // namespace ant::reconstruct
